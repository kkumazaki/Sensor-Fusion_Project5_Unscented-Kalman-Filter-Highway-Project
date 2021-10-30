# SFND_Unscented_Kalman_Filter
Sensor Fusion UKF Highway Project Starter Code

<img src="media/ukf_highway_tracked.gif" width="700" height="400" />

## Goal
In this project I implemented an Unscented Kalman Filter to estimate the state of multiple cars on a highway using noisy lidar and radar measurements. Passing the project requires obtaining RMSE values that are lower that the tolerance outlined in the project rubric as below:

- px, py, vx, vy output coordinates must have an   RMSE <= [0.30, 0.16, 0.95, 0.70] after running for longer than 1 second.

## Overview
### (1)Source Codes
- main.cpp: main code of this project.  
  I modified the argument of "stepHighway()" to log the RMSE.
- measurement_package.h: Contains the class MeasurementPackage.
- tools.h, tools.cpp: Contains useful functions such as noise(), lidarSense(), radarSense(), ukfResults(), CalculateRMSE(), savePcd(), loadPcd().
- highway.h: Handle logics for creating traffic on highway and animating it. We can select the settings as following:  
  * trackCars: choose cars to measure. max: 3.
  * visualize_lidar: show red sphere
  * visualize_radar: show pink lines
  * visualize_pcd: show gray point clouds
  * projectTime, projectSteps: show green spheres that predict the vehicle paths  

- ukf.h, ukf.cpp: main project codes that I modified by reffering to the Lesson codes.


### (2)Points
- I modified process noise standard deviations.  
  Default values: 30 were too large, so I modified them to the moderate values.
  std_a_ = 1.0;  
  std_yawdd_ = 0.3*M_PI;
- I set the initial state vector = 0.
- At first I set the initial covariance matrix = identity matrix, but the initial error of velocity was large and overshooted the target RMSE. So I changed the covariance of v and phi as following:  
  P_(2,2) = 100;  
  P_(3,3) = 100;  

## Result
### (1)Sensor fusion of Lidar and Radar
The following is the resulting RMSE of Sensor fusion of Lidar and Radar.
It achieved the target goal of project rubric.
<img src="result/result_Lidar_Radar.png" />

### (2)Only Lidar
The following is the resultig RMSE of only Lidar measurement.
It didn't achieve the target goal of project rubric. The accuracy of position(x and y) look fine, but the accuracy of velocity(vx and vy) are not good.
<img src="result/result_Lidar.png" />

### (3)Only Radar
The following is the resultig RMSE of only Radar measurement.
It didn't achieve the target goal of project rubric. The accuracy of both position and velocity are not good.
<img src="result/result_Radar.png" />


## SFND_Unscented_Kalman_Filter program
The main program can be built and ran by doing the following from the project top directory.

1. mkdir build
2. cd build
3. cmake ..
4. make
5. ./ukf_highway

Note that the programs that need to be written to accomplish the project are src/ukf.cpp, and src/ukf.h

The program main.cpp has already been filled out, but feel free to modify it.

<img src="media/ukf_highway.png" width="700" height="400" />

`main.cpp` is using `highway.h` to create a straight 3 lane highway environment with 3 traffic cars and the main ego car at the center. 
The viewer scene is centered around the ego car and the coordinate system is relative to the ego car as well. The ego car is green while the 
other traffic cars are blue. The traffic cars will be accelerating and altering their steering to change lanes. Each of the traffic car's has
it's own UKF object generated for it, and will update each indidual one during every time step. 

The red spheres above cars represent the (x,y) lidar detection and the purple lines show the radar measurements with the velocity magnitude along the detected angle. The Z axis is not taken into account for tracking, so you are only tracking along the X/Y axis.

---

## Other Important Dependencies
* cmake >= 3.5
  * All OSes: [click here for installation instructions](https://cmake.org/install/)
* make >= 4.1 (Linux, Mac), 3.81 (Windows)
  * Linux: make is installed by default on most Linux distros
  * Mac: [install Xcode command line tools to get make](https://developer.apple.com/xcode/features/)
  * Windows: [Click here for installation instructions](http://gnuwin32.sourceforge.net/packages/make.htm)
* gcc/g++ >= 5.4
  * Linux: gcc / g++ is installed by default on most Linux distros
  * Mac: same deal as make - [install Xcode command line tools](https://developer.apple.com/xcode/features/)
  * Windows: recommend using [MinGW](http://www.mingw.org/)
 * PCL 1.2

## Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./ukf_highway`

## Editor Settings

We've purposefully kept editor configuration files out of this repo in order to
keep it as simple and environment agnostic as possible. However, we recommend
using the following settings:

* indent using spaces
* set tab width to 2 spaces (keeps the matrices in source code aligned)

## Code Style

Please stick to [Google's C++ style guide](https://google.github.io/styleguide/cppguide.html) as much as possible.

## Generating Additional Data

This is optional!

If you'd like to generate your own radar and lidar modify the code in `highway.h` to alter the cars. Also check out `tools.cpp` to
change how measurements are taken, for instance lidar markers could be the (x,y) center of bounding boxes by scanning the PCD environment
and performing clustering. This is similar to what was done in Sensor Fusion Lidar Obstacle Detection.

## Project Instructions and Rubric

This information is only accessible by people who are already enrolled in Sensor Fusion. 
If you are enrolled, see the project page in the classroom
for instructions and the project rubric.
