f = read_artsvar('example1','f_mono');
z = read_artsvar('example1','z_abs');  
abs_tg = read_artsvar('example1','abs_per_tg');  
mesh(z,f,log(abs_tg{1}))    
