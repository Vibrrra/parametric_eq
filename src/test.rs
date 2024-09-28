


use std::{f32::INFINITY, rc::Rc};

use biquad::{Coefficients, Hertz};
use slint::{ComponentHandle, Model, ModelRc, VecModel, Weak};
use crate::{signal_analysis::{calculate_biquad_response, f_to_lin_skew, points_to_smooth_svg_path}, ui::*, EqCurves}; //{EQManagerUI,TestWindow, EqFilter, EQGraphManager};


fn map_filtertype(filter_type: EQFilterTypeEnum, gain: f32) -> biquad::Type<f32> {
    match filter_type {
        
        EQFilterTypeEnum::Lowpass => biquad::Type::LowPass,
        EQFilterTypeEnum::Highpass => biquad::Type::HighPass,
        EQFilterTypeEnum::Bandpass => biquad::Type::BandPass,
        EQFilterTypeEnum::Peaking=> biquad::Type::PeakingEQ(gain),
        EQFilterTypeEnum::Lowshelf => biquad::Type::LowShelf(gain),
        EQFilterTypeEnum::Highshelf => biquad::Type::HighShelf(gain),
        EQFilterTypeEnum::Notch => biquad::Type::Notch,
    }
}


pub struct AppHandler {
    handler: Option<TestWindow>,
    eq_control: EqControl,
} 

impl AppHandler {
    
    pub fn new() -> Self {      
        AppHandler {
            handler: None,
            eq_control: EqControl {  },
        
        }
    }   

    pub fn init(&mut self){
        let window = TestWindow::new().expect("Cannot create main window!");
        let (tx, rx) = std::sync::mpsc::channel();
        self.eq_control.init(&window, tx);
        self.handler = Some(window);
    }

    pub fn run(&mut self) {
        let window = self.handler.as_ref().unwrap();
        window.run();
    }    
}  


// EQManagerSettings.eq_graph_height/2px - (EQManagerSettings.eq_graph_height/2px * (float / EQManagerSettings.max_gain));
pub fn gain_to_my(gain: f32, eq_graph_height: f32, max_gain: f32) -> f32 {
    eq_graph_height/2.0 - (eq_graph_height/2.0 * (gain / max_gain))

}




pub struct EqControl {}

impl EqControl {
    
    pub fn init(&self, handle: &TestWindow, tx: std::sync::mpsc::Sender<Coefficients<f32>>) {
        let eq_manager = handle.global::<EQManagerUI>();
        let eq_settings = handle.global::<EQManagerSettings>();
        
        let mut eq_curves = EqCurves {
            curves: vec![0.0f32; eq_settings.get_eq_filters().row_count() * eq_settings.get_num_freq_points() as usize],
            sum_curve: vec![0.0; eq_settings.get_num_freq_points() as usize],
            n_curves: eq_settings.get_eq_filters().row_count(),
            n_samples: eq_settings.get_num_freq_points() as usize,
            is_active: [true; 8],
        };
        
        
        eq_manager.set_test(32.0);

        eq_manager.on_init_eq_filters({
            let ww = handle.as_weak();

            move || {
                let wu = ww.upgrade().unwrap();
                
                let eq_filter: Vec<EqFilter> = Vec::new();
                let eq_filter = Rc::new(VecModel::from(eq_filter));      
                let n_filters = eq_filter.row_count();
                let init_filters = wu.global::<EQManagerSettings>().get_eq_filters();
                init_filters.iter()
                    .enumerate()
                    .for_each(|(n, mut init_filter)| {
                        
                        // scale frequency to canvas coords
                        let x = f_to_lin_skew(
                            wu.global::<EQManagerSettings>().get_min_freq(), 
                            wu.global::<EQManagerSettings>().get_max_freq(), 
                            wu.global::<EQManagerSettings>().get_num_freq_points() as usize, 
                            wu.global::<EQManagerSettings>().get_skew_factor(), init_filter.frequency) ;
                        init_filter.circle_x = wu.global::<EQManagerSettings>().get_eq_graph_width() 
                            * x / wu.global::<EQManagerSettings>().get_num_freq_points() as f32;

                        // scale gain to canvas coords
                        let y = gain_to_my(
                            init_filter.gain,
                            wu.global::<EQManagerSettings>().get_eq_graph_height(),
                            wu.global::<EQManagerSettings>().get_max_gain() );
                        init_filter.circle_y = y;    
                        
                        eq_filter.push(init_filter.clone());    
                        
                    }) ;

                let eq_filter = ModelRc::new((eq_filter.clone()));      
                
                wu.global::<EQManagerUI>().set_eq_filters(eq_filter);    
                for k in 0..n_filters {
                    wu.global::<EQManagerUI>().invoke_calculate_filter_coefficients(k as i32);
                }
                // wu.global::<EQManagerUI>().invoke_calculate_filter_coefficients(0);
            }
        });
        
        eq_manager.on_set_filter_frequency({
            let ww = handle.as_weak();
            move |id, freq| { 
                let wu = ww.upgrade().unwrap();
                let filters  = wu.global::<EQManagerUI>().get_eq_filters();
                let mut filter = filters.row_data(id as usize).unwrap();
                
                filter.frequency = freq.clamp(wu.global::<EQManagerSettings>().get_min_freq(), wu.global::<EQManagerSettings>().get_max_freq());
                filter.huge_frequency = freq > 10000.0;
                filters.set_row_data(id as usize, filter);
                
        }}
        );

        eq_manager.on_set_filter_gain({
            let ww = handle.as_weak();
            move |id, gain| { 
                let wu = ww.upgrade().unwrap();
                let filters  = wu.global::<EQManagerUI>().get_eq_filters();
                let mut filter = filters.row_data(id as usize).unwrap();
                filter.gain = gain.clamp(-wu.global::<EQManagerSettings>().get_max_gain(), wu.global::<EQManagerSettings>().get_max_gain());
                filters.set_row_data(id as usize, filter);
        }}
        );
        eq_manager.on_set_filter_q({
            let ww = handle.as_weak();
            move |id, q| { 
                let wu = ww.upgrade().unwrap();
                let filters  = wu.global::<EQManagerUI>().get_eq_filters();
                let mut filter = filters.row_data(id as usize).unwrap();
                filter.q = q.clamp(wu.global::<EQManagerSettings>().get_min_q(), wu.global::<EQManagerSettings>().get_max_q());
                filters.set_row_data(id as usize, filter);
                
        }}
        );
        eq_manager.on_set_filter_type({
            let ww = handle.as_weak();
            move |id, filter_type| { 
                let wu = ww.upgrade().unwrap();
                let filters  = wu.global::<EQManagerUI>().get_eq_filters();
                let mut filter = filters.row_data(id as usize).unwrap();
                filter.filter_type = filter_type;
                filters.set_row_data(id as usize, filter);
                wu.global::<EQManagerUI>().invoke_calculate_filter_coefficients(id);
        }}
        );
        eq_manager.on_set_filter_enabled({
            let ww = handle.as_weak();
            move |id, enabled| { 
                let wu = ww.upgrade().unwrap();
                let filters  = wu.global::<EQManagerUI>().get_eq_filters();
                let mut filter = filters.row_data(id as usize).unwrap();
                filter.enabled = enabled;
                filters.set_row_data(id as usize, filter);
        }}
        );

        // nearest neighbor
        eq_manager.on_select_nearest_draggable({
            let ww = handle.as_weak();
            move|x,y| {
                let ui = ww.upgrade().unwrap();
                let filters = ui.global::<EQManagerUI>().get_eq_filters();
        
                let min_val: f32 = INFINITY;
                let mut sel: i32 = -1;
                let min_val: f32 = filters.iter().enumerate().fold(INFINITY, |min, (n, filter)| { 
                    let dist: f32 = ((filter.circle_x - x).powi(2) + (filter.circle_y - y).powi(2)).sqrt();
                    if dist < min && dist < 100.0{
                        sel = n as i32; dist
                    } else {
                        min
                    }
        
                });
                // if sel != -1 {
                    ui.global::<EQManagerUI>().set_selected_filter(sel as i32);
                // }
                
           
            }
        });

        eq_manager.on_calculate_filter_coefficients({
            let ww = handle.as_weak();

            move |id| {
                let wu =ww.upgrade().unwrap();
                let fs = wu.global::<EQManagerSettings>().get_sample_rate();
                let filter = wu.global::<EQManagerUI>().get_eq_filters().row_data(id as usize).unwrap(); 
                let filter_type_biquad = map_filtertype(filter.filter_type.r#type, filter.gain);
                let gain = filter.gain;
                let freq = filter.frequency;
                let q = filter.q;
                let filter_type = filter.filter_type;//map_filtertype(filter.filter_type, gain);
                let coeffs = Coefficients::<f32>::from_params(
                    filter_type_biquad, 
                    Hertz::<f32>::from_hz(fs).unwrap(), 
                    Hertz::<f32>::from_hz(freq).unwrap(),
                    q
                ).unwrap();

                let (filter_response, _) = 
                    calculate_biquad_response(
                        (coeffs.b0,coeffs.b1,coeffs.b2), 
                        (1.0, coeffs.a1, coeffs.a2), 
                        fs,
                        wu.global::<EQManagerSettings>().get_num_freq_points() as usize);
                

                // update sum curve
                eq_curves.copy_from_slice(id as usize, &filter_response, filter.enabled);
                eq_curves.update();
                let points_sum_curve: Vec<(f32,f32)> = 
                eq_curves.sum_curve
                    .iter()
                    .enumerate()
                    .map(|x| {
                        // map_to_canvas(x, wu.global::<EQManagerSettings>()., n_freq_points, eq_canvas_half_height)
                        (x.0 as f32 
                            / wu.global::<EQManagerSettings>().get_num_freq_points() as f32 
                            * wu.global::<EQManagerSettings>().get_eq_graph_width(),
                            /* wu.global::<EQManagerSettings>().get_max_gain() */ - x.1 
                            / wu.global::<EQManagerSettings>().get_max_gain() 
                            * wu.global::<EQManagerSettings>().get_eq_graph_height() / 2.0 )
                    }).collect();        
                let svg_sum = points_to_smooth_svg_path(&points_sum_curve, 0.2);


                let points: Vec<(f32,f32)> = 
                    filter_response
                        .iter()
                        .enumerate()
                        .map(|x| {
                            // map_to_canvas(x, wu.global::<EQManagerSettings>()., n_freq_points, eq_canvas_half_height)
                            (x.0 as f32 
                                / wu.global::<EQManagerSettings>().get_num_freq_points() as f32 
                                * wu.global::<EQManagerSettings>().get_eq_graph_width(),
                                /* wu.global::<EQManagerSettings>().get_max_gain() */ - x.1 
                                / wu.global::<EQManagerSettings>().get_max_gain() 
                                * wu.global::<EQManagerSettings>().get_eq_graph_height() / 2.0 )
                        }).collect();

                let svg = points_to_smooth_svg_path(&points, 0.2);
                
                let mut filter = wu.global::<EQManagerUI>().get_eq_filters().row_data(id as usize).unwrap();
                filter.curve = svg;
                wu.global::<EQManagerUI>().get_eq_filters().set_row_data(id as usize, filter);
                wu.global::<EQManagerUI>().set_eq_sum_curve(svg_sum);
                // send to audiothread         
                let _ = tx.send(coeffs);                        
            }
        });

        eq_manager.invoke_init_eq_filters();
    }
}

fn do_smthing(handle: &Weak<TestWindow>) {}