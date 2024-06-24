#include <memory>
#include <chrono>
#include <math.h>
#include <iostream>

using namespace std::chrono;

inline static constexpr double k_d = 2.0; // net force reaction parameter
inline static constexpr double k_q = 2.0; // torque reaction parameter
inline static constexpr double g = 9.81; // convert from m/s^2 to m/ns^2
inline static constexpr double max_thrust = 10.0;

struct drone_state
{
    struct angle_t 
    {
        double theta, phi, alpha; // angle about z, x, y
    };

    struct vector_t
    {
        double x, y, z;
    };

    struct rod
    {
        double Fr, Fl;
    } rod_a, rod_b;

    angle_t  angle, omega; // omega in rad/ns
    vector_t position, velocity, target_pos; // in m and m/ns
    vector_t d, dddt; 

    drone_state(double m, double _I, double _r) :
        mass(m), I(_I), r(_r)
    {   }
    
    const double mass, r, I;
};

double clip(double val, double min, double max)
{
    return ( val < max ? ( val > min ? val : min ) : max );
}

double torque(
    const drone_state& state, 
    double angle, 
    double angle_dot, 
    double Fr, 
    double Fl, 
    double Fd_x)
{
    const auto theta_parallel = ( Fr + Fl != 0 ? asin( clip(Fd_x, -( Fr + Fl ), Fr + Fl) / (Fr + Fl) ) : 0.0 );
    return -2.0 * sqrt(state.I * k_q) * angle_dot - k_q * ( angle - theta_parallel );
}

int main()
{
    drone_state state(1.0, 0.75, 2.0);
    state.target_pos.y = 1.0;

    //adc_gpio_init( analog_pin )

    auto start = high_resolution_clock::now();
    while (true)
    {
        auto dt = duration<double, seconds::period> (start - high_resolution_clock::now()).count();
        start = high_resolution_clock::now();

        // Gryoscope
        {
            // read in gyro data as deg/s
            // convert to rad/ns

            //state.omega.theta = input
            //state.omega.phi = input
            //state.omega.alpha = input

            state.angle.theta += state.omega.theta * dt;
            state.angle.phi   += state.omega.phi   * dt;
            state.angle.alpha += state.omega.alpha * dt;
        }

        // Precompute the trigonometric functions
        const double cos_theta = cos(state.angle.theta), sin_theta = sin(state.angle.theta);
        const double cos_phi   = cos(state.angle.phi  ), sin_phi   = sin(state.angle.phi  );
        const double cos_alpha = cos(state.angle.alpha), sin_alpha = sin(state.angle.alpha);

        // Accelerometer
        {
            // read in accelerometer data (probably in units of g) convert to m/ns^2
            double ap_x, ap_y, ap_z; // m/s^2

            // perform a coordinate transformation from drone space to world space
            const double a_x = ap_x * ( sin_alpha * sin_theta * sin_phi + cos_alpha * cos_theta ) + 
                         ap_y * sin_theta * cos_phi + 
                         ap_z * ( cos_alpha * sin_theta * sin_phi - sin_alpha * cos_theta );
            
            const double a_y = ap_x * ( sin_alpha * cos_theta * sin_phi - cos_alpha * sin_theta ) +
                               ap_y * cos_theta * cos_phi +
                               ap_z * ( cos_alpha * cos_theta * sin_phi + sin_alpha * sin_theta );

            const double a_z = ap_x * sin_alpha * cos_phi -
                               ap_y * sin_phi +
                               ap_z * cos_alpha * cos_phi;

            state.velocity.x += a_x * dt; state.position.x += state.velocity.x * dt;
            state.velocity.y += a_y * dt; state.position.x += state.velocity.y * dt;
            state.velocity.z += a_z * dt; state.position.x += state.velocity.z * dt;
        }

        // Using this information find the displacement vector 
        // We could possibly just calculate the displacement velocity 
        // by projecting the world_space velocity onto the displacement unit vector
        {
            const auto d = drone_state::vector_t{
                state.target_pos.x - state.position.x,
                state.target_pos.y - state.position.y,
                state.target_pos.z - state.position.z
            };

            state.dddt = drone_state::vector_t{
                (d.x - state.d.x) / dt,
                (d.y - state.d.y) / dt,
                (d.z - state.d.z) / dt
            };
            state.d = d;
        }

        // Calculate the desired force
        const auto coeff_dddt = 2.0 * sqrt(state.mass * k_d);
        const auto Fd = drone_state::vector_t{
            coeff_dddt * state.dddt.x + k_d * state.d.x,
            coeff_dddt * state.dddt.y + k_d * state.d.y,
            coeff_dddt * state.dddt.z + k_d * state.d.z
        };

        // rod a is in the xy-drone plane
        // rod b is in the zy-drone plane
        // We need to take our target force and split it into the plane aligned
        // with rod a and the plane aligned with rod b. Then we run our computations
        // to generate torques for each and therefore the individual thrusts for all 
        // 4 propellors

        // basis vectors for a and b start by rotating world basis about y-axis by alpha degrees
        // then, for a, we rotate about the x-axis by phi degrees
        // we use this transformation to find the components of Fd in this plane (ignore z?)
        // then we compute torque for rod a
        {
            const auto Fd_ax = Fd.x * cos_alpha - Fd.z * sin_alpha;
            const auto Fd_ay = Fd.x * sin_alpha * sin_phi + Fd.y * cos_phi + Fd.z * cos_alpha * sin_phi;
            //const auto Fd_az = Fd.x * sin_alpha * cos_phi - Fd.y * sin_phi + Fd.z * cos_alpha * cos_phi;

            const auto T_rod_a = torque(
                state, 
                state.angle.theta, 
                state.omega.theta, 
                state.rod_a.Fr, 
                state.rod_a.Fl, 
                Fd_ax);

            const auto gravity_force_thrust = (state.mass * g + Fd_ay) / ( 2.0 * cos_theta );
            state.rod_a.Fr = clip( gravity_force_thrust + T_rod_a / ( 2.0 * state.r ), 0, max_thrust );
            state.rod_a.Fl = clip( gravity_force_thrust - T_rod_a / ( 2.0 * state.r ), 0, max_thrust );
        }

        // we take the basis vectors, but this time rotate about the z axis by theta degrees
        // then we use this transformation to find the components of Fd in this plane (ignore x?)
        // then we compute torque for rod b
        {
            //const auto Fd_ax = Fd.x * cos_alpha * cos_theta + Fd.y * sin_theta - Fd.z * cos_alpha * cos_theta;
            const auto Fd_ay = -1.0 * Fd.x * cos_alpha * sin_theta + Fd.y * cos_theta + Fd.z * sin_alpha * sin_theta;
            const auto Fd_az = Fd.x * sin_alpha + Fd.z * cos_alpha;

            const auto T_rod_a = torque(
                state, 
                state.angle.phi, 
                state.omega.phi, 
                state.rod_b.Fr, 
                state.rod_b.Fl, 
                Fd_az);

            const auto gravity_force_thrust = (state.mass * g + Fd_ay) / ( 2.0 * cos_phi );
            state.rod_b.Fr = clip( gravity_force_thrust + T_rod_a / ( 2.0 * state.r ), 0, max_thrust );
            state.rod_b.Fl = clip( gravity_force_thrust - T_rod_a / ( 2.0 * state.r ), 0, max_thrust );
        }

        // Note: after the transformations, we still have a 3d vector, but we are working
        //       in a plane... so do we just ignore the (x, z)-dimension?

        // At this point each engine has a thrust vector, so it's just a matter of converting 
        // these values to voltages and writing them out to each servo

        //adc_select_input( analog_pin )
        //uint16_t val = adc_read()

        // output PWM
        // You set the cycle period and then tell a PIN to write high for N cycles
        // (cycles per period) / processor clock cycle = time of period in seconds
        //https://github.com/raspberrypi/pico-examples/blob/master/pwm/hello_pwm/hello_pwm.c
    }

}