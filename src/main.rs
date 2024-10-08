
// #![windows_subsystem = "windows"]

use std::{f32::{consts::E, INFINITY, NAN, NEG_INFINITY}, thread, vec};

use biquad::{frequency, Coefficients, Hertz, Type};
use filters::BiquadFilter;
use signal_analysis::{calculate_biquad_response, create_grid_svg, f_to_lin_skew, lin_to_f, lin_to_f_skew, points_to_smooth_svg_path};
use slint::{ComponentHandle, Model, ModelExt, ModelRc, SharedString, VecModel};

pub mod UI{ slint::include_modules!();}
mod filters;
mod signal_analysis;
// #[derive(Clone)]
// struct EQCurveData {
//     svg: SharedString,
//     highlighted: bool,
//     control_points: (f32, f32),
// }

#[derive(Clone)]
pub struct EqCurves {
    curves: Vec<f32>,
    n_samples: usize,
    n_curves: usize,
    is_active: [bool; N_FILTERS],
    sum_curve: Vec<f32>
}
impl EqCurves {
    pub fn copy_from_slice(&mut self, n: usize, slice: &[f32], active: bool) {
        assert!(slice.len() == self.n_samples);
        assert!(n < self.n_curves);
        self.is_active[n] = active;
        let o = &mut self.curves[n*self.n_samples .. (n+1)*self.n_samples];
        o.copy_from_slice(slice);
    }

    // sum all active eq_curves
    pub fn update(&mut self) {
        // flush sum
        self.sum_curve.iter_mut().for_each(|x| *x = 0.0);
        for i in 0..self.n_curves {            
            if self.is_active[i] {
                // sum the curves
                self.sum_curve.iter_mut()
                .zip(self.curves[i*self.n_samples .. (i+1)*self.n_samples].iter())
                .for_each(|(a, b)| *a += b);
                }
        }
    }
}

fn match_biquad_type(biquad_type: UI::FilterType, gain: f32) -> Type<f32> {
    match biquad_type {
        UI::FilterType::LowPass => Type::LowPass,
        UI::FilterType::HighPass => Type::HighPass,
        UI::FilterType::BandPass => Type::BandPass,
        UI::FilterType::Notch => Type::Notch,   
        UI::FilterType::PeakingEQ => Type::PeakingEQ(gain),
        UI::FilterType::LowShelf => Type::LowShelf(gain),
        UI::FilterType::HighShelf => Type::HighShelf(gain),
    }   
}


// build svg grid

const N_FILTERS: usize = 8;
fn main() {
    println!("Hello, world!");

    let ui = UI::AppWindow::new().unwrap();
    let ui_weak = ui.as_weak().unwrap();

    let s = "";
    let s = SharedString::from(s);
    let highlighted = false;
    let cp = (0.0f32, 0.0f32);
    
    let sampling_rate = 48000.0;
    let sampling_rate_hz = Hertz::<f32>::from_hz(sampling_rate).unwrap();
    let n_freq_points = 250;
    
    let mut eq_curves = EqCurves {
        curves: vec![0.0f32; n_freq_points*N_FILTERS],
        sum_curve: vec![0.0; n_freq_points],
        n_curves: N_FILTERS,
        n_samples: n_freq_points,
        is_active: [true; N_FILTERS],
    };
    // let (tx, rx) = crossbeam::channel::bounded(1);
    

    // Function that activates a based on the Mouse-Click-Position on the EQ-Canvas
    //  - clickable_radius is the radius around the mouse click that is considered as a click on a filter
    let clickable_radius = 10.0; 
    ui.global::<UI::DraggableLogic>().on_get_nearest_filter(move|d| {
        let n = d.row_count();
        // println!("N: {n}");
        let min_val = d.iter().fold(INFINITY, |min, val| { if val.dist < min {val.dist} else {min}});
        // println!("Min: {min_val}");
        let index = d.iter().position(|x| x.dist == min_val).unwrap();
        // dbg!(index as i32);
        if min_val < clickable_radius {
            index as i32
        } else {
            -1
        }
        // min_index
    });


    // Pre calculate some values 
    // TODO: Delete this and make it dynamic in the future
    let eq_canvas_half_width = ui.global::<UI::EQGraphManager>().get_eq_graph_width() / 2.0;
    let eq_canvas_half_height = ui.global::<UI::EQGraphManager>().get_eq_graph_height() / 2.0;
    let n_freq_points = ui.global::<UI::EQGraphManager>().get_n_freq_points();
    let skew_f = ui.global::<UI::EQGraphManager>().get_skew_factor();
    
    let ui_c = ui.clone_strong();
    ui.global::<UI::EqCanvasGrid>().on_create_grid( move||{
            let fp = ui_c.global::<UI::EqCanvasGrid>().get_freq_points();
            let gp: ModelRc<f32> = ui_c.global::<UI::EqCanvasGrid>().get_gain_points();
            let mut h = 
                fp.iter()
                    .map(|x| {x
                        // f_to_lin_skew(20.0, 20000.0, n_freq_points as usize, skew_f, x)
                    })
                    .collect::<Vec<(f32)>>();
            let mut g = 
                gp.iter()
                    .map(|x| {x
                        // f_to_lin_skew(20.0, 20000.0, n_freq_points as usize, skew_f, x)
                    })
                    .collect::<Vec<(f32)>>();
            // h.iter_mut().for_each(|x| *x = map_to_canvas((0, x), eq_canvas_half_width, n_freq_points, eq_canvas_half_height).1);
            let svg_f: Vec<(f32)> =  h.iter().enumerate().map(|(u, x)| {
                let temp = 400.0 * f_to_lin_skew(20.0, 20000.0, n_freq_points as usize, skew_f, (*x)) / n_freq_points as f32; // 250.0 before
                // let o = map_to_canvas(( temp,x.1), eq_canvas_half_width, n_freq_points, eq_canvas_half_height)
                temp
            }).collect();
            let svg_g: Vec<f32> = g.iter().map(|x| {
                // eq_canvas_half_height - (x.1 * eq_canvas_half_height / 30.0)}).collect();
                eq_canvas_half_height - eq_canvas_half_height * ((x +  30.0)/60.0) }).collect();
            let grid = create_grid_svg(svg_f, svg_g);//svg_g);
            
            ui_c.global::<UI::EqCanvasGrid>().set_svg(grid);
        }
    );

    let ui_c = ui.clone_strong();
    ui.global::<UI::EQCanvasLogic>().on_calc_new_filter(move |filter| {
        
        // get current scaling factor
        // let db_scaling_factor = ui_c.global::<UI::EQGraphManager>().get_db_scaling_factor();
        
        // scale gain. Mouse-Y-Pos -> Gain
        let gain = 30.0  * (eq_canvas_half_height - filter.y) / eq_canvas_half_height;//(2.0 * db_scaling_factor);
        // dbg!(gain);
        // scale frequency. Mouse-X-Pos -> Frequency
        let freq = lin_to_f_skew(20.0, 20000.0, eq_canvas_half_width as usize * 2 , skew_f, filter.x);
        
        // get filter type. Converts UI-Enum-Variant(Slint) to Biquad-Enum-Variant (Rust)
        // let fl: UI::FilterType = filter.filter_type;
        let filter_type = match_biquad_type(filter.filter_type.1, gain); 
        
        let q_value = filter.q_value.clamp( 
            ui_c.global::<UI::EQGraphManager>().get_min_q(),
            ui_c.global::<UI::EQGraphManager>().get_max_q());
        
        // calculate new coefficients
        let new_coeffs = Coefficients::<f32>::from_params(filter_type, sampling_rate_hz, Hertz::<f32>::from_hz(freq).unwrap(), q_value).unwrap();
        
        // calculate new response
        let (new_response, _) = signal_analysis::calculate_biquad_response(
            (new_coeffs.b0, new_coeffs.b1, new_coeffs.b2), 
            (1.0, new_coeffs.a1, new_coeffs.a2), 
             sampling_rate, n_freq_points as usize);

        // update eq_curves container
        eq_curves.copy_from_slice(filter.id as usize, &new_response, filter.active);
        
        // update sum curve
        eq_curves.update();

        // update eq_graph_manager with new sum curve
        // let eq_cruves_vector = &eq_curves.sum_curve.iter().enumerate().map(|x| {(x.0 as f32 *eq_canvas_half_width / (n_freq_points-1) as f32, 50.0 - *x.1 * db_scaling_factor )}).collect::<Vec<(f32, f32)>>();
        // fn map_to_canvas(x: (usize, &f32), eq_canvas_half_width: f32, n_freq_points: i32, eq_canvas_half_height: f32) -> (f32, f32) {
        let eq_cruves_vector = &eq_curves.sum_curve.iter().enumerate().map(|x| {map_to_canvas(x, eq_canvas_half_width, n_freq_points, eq_canvas_half_height)}).collect::<Vec<(f32, f32)>>();
        let smooth_svg_path = points_to_smooth_svg_path(&eq_cruves_vector, 0.5);
        ui_c.global::<UI::EQGraphManager>().set_eq_sum_curve(smooth_svg_path);

        // update eq_canvas with new curve
        // let biquad_response_vector: Vec<(f32, f32)> = new_response.iter().enumerate().map(|x| {(x.0 as f32 *eq_canvas_half_width / (n_freq_points-1) as f32, 50.0 - *x.1 * db_scaling_factor )}).collect();
        // let biquad_response_vector: Vec<(f32, f32)> = new_response.iter().enumerate().map(|x| {(x.0 as f32 *eq_canvas_half_width / (n_freq_points-1) as f32, eq_canvas_half_height/2.0 - (*x.1 /30.0) * (eq_canvas_half_height/2.0) )}).collect();
        let biquad_response_vector: Vec<(f32, f32)> = new_response.iter().enumerate().map(|x| {map_to_canvas(x, eq_canvas_half_width, n_freq_points, eq_canvas_half_height)}).collect();
        let smooth_svg_path = points_to_smooth_svg_path(&biquad_response_vector, 0.5);
        ui_c.global::<UI::EQGraphManager>().set_new_curve(smooth_svg_path); 
        // return new curve (svg path)
        // smooth_svg_path
    });
    
   
    let ui_clone = ui.clone_strong();

    let sampling_rate_hz = Hertz::<f32>::from_hz(48000.0).unwrap();
    let t = thread::spawn(move|| {
        
        let mut biquads =vec![BiquadFilter{
            biquad_coeffs: Coefficients::<f32>::from_params(Type::PeakingEQ(3.0), sampling_rate_hz, Hertz::<f32>::from_hz(250.0).unwrap(), 2.0).unwrap(),
            biquad_type: Type::<f32>::LowPass,
            gain: 0.0,
            freq: 0.0,
        }];

        loop {
            let coeffs = &mut biquads[0].biquad_coeffs;
            let new_coeffs = Coefficients::<f32>::from_params(Type::PeakingEQ(3.0), sampling_rate_hz, Hertz::<f32>::from_hz(250.0).unwrap(), 2.0).unwrap();
            *coeffs = new_coeffs;
            thread::sleep(std::time::Duration::from_millis(10));
        }

    });
    ui.run();

}

fn map_to_canvas(x: (usize, &f32), eq_canvas_half_width: f32, n_freq_points: i32, eq_canvas_half_height: f32) -> (f32, f32) {
    (x.0 as f32 * eq_canvas_half_width / (n_freq_points-1) as f32, eq_canvas_half_height/2.0 - (*x.1 /30.0) * (eq_canvas_half_height/2.0) )
}


#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_f_to_lin_skew() { 
        let n_freq_points = [50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000];
        let skew_f = 0.3;
        for n in n_freq_points.iter() {
            let h = 800.0 * f_to_lin_skew(20.0, 20000.0, 250, skew_f, *n as f32) / 250.0;
            dbg!(h);
        }
    }
}