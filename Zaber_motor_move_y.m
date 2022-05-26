function Zaber_motor_move_y(motor_device, y_amount, y_units)
motor_tic = tic;
motor_axis = motor_device.getAxis(1);
% axis.home();
% Move to the 10mm position
%axis.moveAbsolute(10, Units.LENGTH_MILLIMETRES);
% Move by an additional 5mm    
motor_axis.moveRelative(y_amount, y_units);    
motor_toc = toc(motor_tic);
disp(['Finished moving motor with elapsed time: ', num2str(motor_toc)]);
end
