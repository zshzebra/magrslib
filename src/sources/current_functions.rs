//! Time-varying current functions for dynamic field computation
//!
//! This module provides the foundation for time-dependent current sources,
//! enabling simulation of PWM signals, switching transients, and other
//! dynamic electromagnetic phenomena.

use std::f64::consts::PI;

/// Trait for time-varying current functions I(t)
///
/// Enables dynamic current sources for realistic electromagnetic simulation
/// including PWM signals, switching transients, and AC analysis.
pub trait CurrentFunction: Send + Sync + std::fmt::Debug {
    /// Evaluate current at given time
    ///
    /// # Arguments
    /// * `time` - Time in seconds
    ///
    /// # Returns
    /// Current in Amperes
    fn evaluate(&self, time: f64) -> f64;

    /// Get a descriptive name for this current function
    fn name(&self) -> &'static str;

    /// Clone this current function (for dynamic dispatch)
    fn clone_box(&self) -> Box<dyn CurrentFunction>;
}

// Implement Clone for Box<dyn CurrentFunction>
impl Clone for Box<dyn CurrentFunction> {
    fn clone(&self) -> Self {
        self.clone_box()
    }
}

/// Constant (DC) current source
///
/// Provides steady current for static field analysis.
#[derive(Debug, Clone)]
pub struct ConstantCurrent {
    pub amplitude: f64, // Current in Amperes
}

impl ConstantCurrent {
    pub fn new(amplitude: f64) -> Self {
        Self { amplitude }
    }
}

impl CurrentFunction for ConstantCurrent {
    fn evaluate(&self, _time: f64) -> f64 {
        self.amplitude
    }

    fn name(&self) -> &'static str {
        "Constant"
    }

    fn clone_box(&self) -> Box<dyn CurrentFunction> {
        Box::new(self.clone())
    }
}

/// PWM (Pulse Width Modulated) current source
///
/// Generates square wave current with configurable frequency, duty cycle,
/// and optional rise/fall times for realistic switching behavior.
#[derive(Debug, Clone)]
pub struct PWMCurrent {
    pub frequency: f64,    // Switching frequency in Hz
    pub duty_cycle: f64,   // Duty cycle from 0.0 to 1.0
    pub amplitude: f64,    // Peak current in Amperes
    pub phase: f64,        // Phase offset in radians
    pub rise_time: f64,    // Rise time in seconds (0 = instant)
    pub fall_time: f64,    // Fall time in seconds (0 = instant)
}

impl PWMCurrent {
    /// Create new PWM current with ideal (instantaneous) switching
    pub fn new(frequency: f64, duty_cycle: f64, amplitude: f64) -> Self {
        Self {
            frequency,
            duty_cycle,
            amplitude,
            phase: 0.0,
            rise_time: 0.0,
            fall_time: 0.0,
        }
    }

    /// Create PWM current with finite rise/fall times
    pub fn with_transitions(
        frequency: f64,
        duty_cycle: f64,
        amplitude: f64,
        rise_time: f64,
        fall_time: f64,
    ) -> Self {
        Self {
            frequency,
            duty_cycle,
            amplitude,
            phase: 0.0,
            rise_time,
            fall_time,
        }
    }

    /// Set phase offset in radians
    pub fn with_phase(mut self, phase: f64) -> Self {
        self.phase = phase;
        self
    }
}

impl CurrentFunction for PWMCurrent {
    fn evaluate(&self, time: f64) -> f64 {
        if self.frequency <= 0.0 {
            return 0.0;
        }

        let period = 1.0 / self.frequency;
        let phase_time = (time + self.phase / (2.0 * PI * self.frequency)) % period;
        let normalized_time = phase_time / period;

        let on_duration = self.duty_cycle;

        if self.rise_time <= 0.0 && self.fall_time <= 0.0 {
            // Ideal switching (instant transitions)
            if normalized_time < on_duration {
                self.amplitude
            } else {
                0.0
            }
        } else {
            // Finite rise/fall times
            let rise_fraction = self.rise_time / period;
            let fall_fraction = self.fall_time / period;

            if normalized_time < rise_fraction {
                // Rising edge
                self.amplitude * (normalized_time / rise_fraction)
            } else if normalized_time < on_duration - fall_fraction {
                // High state
                self.amplitude
            } else if normalized_time < on_duration {
                // Falling edge
                let fall_progress = (normalized_time - (on_duration - fall_fraction)) / fall_fraction;
                self.amplitude * (1.0 - fall_progress)
            } else {
                // Low state
                0.0
            }
        }
    }

    fn name(&self) -> &'static str {
        "PWM"
    }

    fn clone_box(&self) -> Box<dyn CurrentFunction> {
        Box::new(self.clone())
    }
}

/// Step current source for transient analysis
///
/// Provides current step at specified time, useful for analyzing
/// switching transients and impulse responses.
#[derive(Debug, Clone)]
pub struct StepCurrent {
    pub step_time: f64,    // Time of step in seconds
    pub amplitude: f64,    // Final current amplitude in Amperes
    pub rise_time: f64,    // Rise time in seconds (0 = instant)
}

impl StepCurrent {
    /// Create instant current step
    pub fn new(step_time: f64, amplitude: f64) -> Self {
        Self {
            step_time,
            amplitude,
            rise_time: 0.0,
        }
    }

    /// Create current step with finite rise time
    pub fn with_rise_time(step_time: f64, amplitude: f64, rise_time: f64) -> Self {
        Self {
            step_time,
            amplitude,
            rise_time,
        }
    }
}

impl CurrentFunction for StepCurrent {
    fn evaluate(&self, time: f64) -> f64 {
        if time < self.step_time {
            0.0
        } else if self.rise_time <= 0.0 {
            // Instant step
            self.amplitude
        } else {
            // Finite rise time
            let elapsed = time - self.step_time;
            if elapsed >= self.rise_time {
                self.amplitude
            } else {
                self.amplitude * (elapsed / self.rise_time)
            }
        }
    }

    fn name(&self) -> &'static str {
        "Step"
    }

    fn clone_box(&self) -> Box<dyn CurrentFunction> {
        Box::new(self.clone())
    }
}

/// Sinusoidal AC current source
///
/// Provides sinusoidal current for AC analysis and frequency domain studies.
#[derive(Debug, Clone)]
pub struct SinusoidalCurrent {
    pub frequency: f64,    // Frequency in Hz
    pub amplitude: f64,    // Peak current in Amperes
    pub phase: f64,        // Phase in radians
    pub dc_offset: f64,    // DC offset in Amperes
}

impl SinusoidalCurrent {
    /// Create sinusoidal current with zero DC offset
    pub fn new(frequency: f64, amplitude: f64) -> Self {
        Self {
            frequency,
            amplitude,
            phase: 0.0,
            dc_offset: 0.0,
        }
    }

    /// Create sinusoidal current with phase offset
    pub fn with_phase(frequency: f64, amplitude: f64, phase: f64) -> Self {
        Self {
            frequency,
            amplitude,
            phase,
            dc_offset: 0.0,
        }
    }

    /// Create sinusoidal current with DC offset (for biased AC analysis)
    pub fn with_dc_offset(frequency: f64, amplitude: f64, dc_offset: f64) -> Self {
        Self {
            frequency,
            amplitude,
            phase: 0.0,
            dc_offset,
        }
    }
}

impl CurrentFunction for SinusoidalCurrent {
    fn evaluate(&self, time: f64) -> f64 {
        self.dc_offset + self.amplitude * (2.0 * PI * self.frequency * time + self.phase).sin()
    }

    fn name(&self) -> &'static str {
        "Sinusoidal"
    }

    fn clone_box(&self) -> Box<dyn CurrentFunction> {
        Box::new(self.clone())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_constant_current() {
        let current = ConstantCurrent::new(5.0);
        assert_relative_eq!(current.evaluate(0.0), 5.0);
        assert_relative_eq!(current.evaluate(100.0), 5.0);
    }

    #[test]
    fn test_pwm_current_ideal() {
        let pwm = PWMCurrent::new(1000.0, 0.5, 10.0); // 1kHz, 50% duty, 10A
        let period = 1.0 / 1000.0; // 1ms period

        // At start of period (high)
        assert_relative_eq!(pwm.evaluate(0.0), 10.0);

        // At 25% of period (still high)
        assert_relative_eq!(pwm.evaluate(0.25 * period), 10.0);

        // At 50% of period (transition to low)
        assert_relative_eq!(pwm.evaluate(0.5 * period), 0.0);

        // At 75% of period (low)
        assert_relative_eq!(pwm.evaluate(0.75 * period), 0.0);

        // At next period start (high again)
        assert_relative_eq!(pwm.evaluate(period), 10.0);
    }

    #[test]
    fn test_pwm_current_duty_cycle() {
        let pwm = PWMCurrent::new(1000.0, 0.25, 10.0); // 25% duty cycle
        let period = 1.0 / 1000.0;

        // High for first 25%
        assert_relative_eq!(pwm.evaluate(0.1 * period), 10.0);
        assert_relative_eq!(pwm.evaluate(0.2 * period), 10.0);

        // Low for remaining 75%
        assert_relative_eq!(pwm.evaluate(0.3 * period), 0.0);
        assert_relative_eq!(pwm.evaluate(0.8 * period), 0.0);
    }

    #[test]
    fn test_step_current() {
        let step = StepCurrent::new(1.0, 5.0); // Step at t=1s to 5A

        assert_relative_eq!(step.evaluate(0.5), 0.0);  // Before step
        assert_relative_eq!(step.evaluate(1.0), 5.0);  // At step
        assert_relative_eq!(step.evaluate(2.0), 5.0);  // After step
    }

    #[test]
    fn test_step_current_with_rise_time() {
        let step = StepCurrent::with_rise_time(1.0, 10.0, 0.1); // 0.1s rise time

        assert_relative_eq!(step.evaluate(0.5), 0.0);   // Before step
        assert_relative_eq!(step.evaluate(1.0), 0.0);   // Start of rise
        assert_relative_eq!(step.evaluate(1.05), 5.0);  // Halfway through rise
        assert_relative_eq!(step.evaluate(1.1), 10.0);  // End of rise
        assert_relative_eq!(step.evaluate(2.0), 10.0);  // After rise
    }

    #[test]
    fn test_sinusoidal_current() {
        let sine = SinusoidalCurrent::new(1.0, 10.0); // 1Hz, 10A amplitude

        assert_relative_eq!(sine.evaluate(0.0), 0.0, epsilon = 1e-10);    // sin(0) = 0
        assert_relative_eq!(sine.evaluate(0.25), 10.0, epsilon = 1e-10);  // sin(π/2) = 1
        assert_relative_eq!(sine.evaluate(0.5), 0.0, epsilon = 1e-10);    // sin(π) = 0
        assert_relative_eq!(sine.evaluate(0.75), -10.0, epsilon = 1e-10); // sin(3π/2) = -1
        assert_relative_eq!(sine.evaluate(1.0), 0.0, epsilon = 1e-10);    // sin(2π) = 0
    }

    #[test]
    fn test_sinusoidal_current_with_dc_offset() {
        let sine = SinusoidalCurrent::with_dc_offset(1.0, 5.0, 2.0); // 2A DC offset

        assert_relative_eq!(sine.evaluate(0.0), 2.0, epsilon = 1e-10);    // 2 + 5*sin(0) = 2
        assert_relative_eq!(sine.evaluate(0.25), 7.0, epsilon = 1e-10);   // 2 + 5*sin(π/2) = 7
        assert_relative_eq!(sine.evaluate(0.75), -3.0, epsilon = 1e-10);  // 2 + 5*sin(3π/2) = -3
    }
}