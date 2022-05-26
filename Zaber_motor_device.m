function [connection, device] = Zaber_motor_device()
% Returns connection and first device. The connection will need to be
% closed manually after usage like so: 
% connection.close();
% if variable is lost, try:
% fclose(instrfindall)
import zaber.motion.Library;
import zaber.motion.ascii.Connection;
import zaber.motion.Units;
Library.enableDeviceDbStore();
connection = Connection.openSerialPort('COM6');
deviceList = connection.detectDevices();
% fprintf('Found %d devices.\n', deviceList.length);
device = deviceList(1);
end

