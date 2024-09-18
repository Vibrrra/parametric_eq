use biquad::Coefficients;

#[derive(Debug, Clone, Copy)]
pub struct BiquadFilter<T> {
    pub biquad_coeffs: Coefficients<f32>,
    pub biquad_type: biquad::Type<T>,
    pub gain: f32,
    pub freq: f32,
}
