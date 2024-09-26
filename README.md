# Implementation of Extended Kalman Filter
The goal is to implement Extended Kalman filter to estimate the pose of a Quadrotor given the IMU and Vicon data. The Vicon being the world frame provides information about the position, orientation, linear and angular velocities. The IMU onboard, through the accelerometer and gyroscope provides the state values. The idea is to make use of the data, use them as control inputs, and determine the state estimates of the quadrotor.
 
