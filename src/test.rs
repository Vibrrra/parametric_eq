


use std::f32::INFINITY;

use biquad::Coefficients;
use slint::{ComponentHandle, Model, Weak};
use crate::ui::*; //{EQManagerUI,TestWindow, EqFilter, EQGraphManager};


pub struct AppHandler {
    handler: Option<TestWindow>,
    eq_control: EqControl,
} 

impl AppHandler {
    pub fn new() -> Self {
        
        // self.weather_display_controller.initialize_ui(&window, self.support_add_city);
        
        
        AppHandler {
            handler: None,
            eq_control: EqControl {  },
        
        }
    }   
    pub fn init(&mut self){
        let window = TestWindow::new().expect("Cannot create main window!");

        self.eq_control.init(&window);

        self.handler = Some(window);
    }

    pub fn run(&mut self) {
        let window = self.handler.as_ref().unwrap();
        window.run();
    }    
}  





pub struct EqControl {
    
}

impl EqControl {
    
    pub fn init(&self, handle: &TestWindow) {
        let eq_manager = handle.global::<EQManagerUI>();
        eq_manager.set_test(32.0);
        eq_manager.on_set_filter_frequency({
            let ww = handle.as_weak();
            move |id, freq| { 
                let wu = ww.upgrade().unwrap();
                let filters  = wu.global::<EQManagerUI>().get_eq_filters();
                let mut filter = filters.row_data(id as usize).unwrap();
                filter.frequency = freq;
                filters.set_row_data(id as usize, filter);
        }}
        );
        eq_manager.on_set_filter_gain({
            let ww = handle.as_weak();
            move |id, gain| { 
                let wu = ww.upgrade().unwrap();
                let filters  = wu.global::<EQManagerUI>().get_eq_filters();
                let mut filter = filters.row_data(id as usize).unwrap();
                filter.gain = gain;
                filters.set_row_data(id as usize, filter);
        }}
        );
        eq_manager.on_set_filter_q({
            let ww = handle.as_weak();
            move |id, q| { 
                let wu = ww.upgrade().unwrap();
                let filters  = wu.global::<EQManagerUI>().get_eq_filters();
                let mut filter = filters.row_data(id as usize).unwrap();
                filter.q = q;
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
        // eq_manager.on_select_nearest_draggable2({
        //     let ww = handle.as_weak();
        //     move |x,y| {
        //         let filters = ww.upgrade().unwrap().global::<EQManagerUI>().get_eq_filters();
        //         let id_selected = 
        //             filters.iter().enumerate().fold(None, |acc, (i, filter)| {
        //                 let x1 = filter.frequency;
        //                 let y1 = filter.gain;
        //                 let x2 = x;
        //                 let y2 = y;
        //                 let dist = (x2 - x1).powi(2) + (y2 - y1).powi(2);
        //                 match acc {
        //                     None => Some((i, dist)),
        //                     Some((_, min_dist)) => {
        //                         if dist < min_dist {
        //                             Some((i, dist))
        //                         } else {
        //                             acc
        //                         }
        //                     }
        //                 }
        //         });
        //     }
        // });
    
    }
}



fn do_smthing(handle: &Weak<TestWindow>) {}