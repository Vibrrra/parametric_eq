


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
        // eq_manager.on_update_freq(|id, freq| {
        //     let o = 1;
        //     let o2 = 2;
        //     let handle_weak = handle.as_weak();
        //     move || {do_smthing(&handle_weak)};
        // });
    
    
    }
}



fn do_smthing(handle: &Weak<TestWindow>) {}