#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"


using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}

int ClosestWaypoint(double x, double y, vector<double> maps_x, vector<double> maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2( (map_y-y),(map_x-x) );

	double angle = abs(theta-heading);

	if(angle > pi()/4)
	{
		closestWaypoint++;
	}

	return closestWaypoint;

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, vector<double> maps_s, vector<double> maps_x, vector<double> maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

// The below four functions are used for building up the next waypoints in the path

vector<double> copyVector(vector<double> v){
    
    vector<double> v2;
    
    for(int i=0;i<v.size();i++){
        v2.push_back(v[i]);
    }
    
    return v2;
}

vector<double> compute_next_d3s(vector<double> path_planning_state){
    double target_velocity = path_planning_state[4];
    double max_acceleration = path_planning_state[5];
    double max_negative_acceleration = max_acceleration * -1;
    double target_jerk = path_planning_state[6];
    
    double current_velocity = path_planning_state[1];
    double current_acceleration = path_planning_state[2];
    
    vector<double> updated_path_planning_state = copyVector(path_planning_state);

    if(abs(current_velocity-target_velocity) < 0.001)
    {
        updated_path_planning_state[3] = 0;
        //return updated_path_planning_state;
    }
    else if(current_velocity < target_velocity)
    {
        if(current_acceleration < max_acceleration)
        {
            updated_path_planning_state[3] = target_jerk;
            //return updated_path_planning_state;
        }
        else
        {
            updated_path_planning_state[3] = 0;
            //return updated_path_planning_state;
        }
    }
    else if (current_velocity > target_velocity){
        if (current_acceleration > max_negative_acceleration)
        {
            updated_path_planning_state[3] = target_jerk * -1;
            //return updated_path_planning_state;
        }
        else
        {
            updated_path_planning_state[3] = 0;
            //return updated_path_planning_state;
        }
    }
    
    
    // TODO: Apply a low pass filter or something similar to dampen the acceleration as the target velocity approaches
    // Something like 'if sufficiently close to the target velocity then negate the jerk causing the acceleration to go back to zero'
    if(abs(current_velocity-target_velocity) < 0.0001){
        updated_path_planning_state[3] = updated_path_planning_state[3] * -1;
    }
    
    
    return updated_path_planning_state;

}

vector<double> compute_next_d2s(vector<double> path_planning_state){
    vector<double> updated_path_planning_state = compute_next_d3s(path_planning_state);

    double next_d2s = path_planning_state[2] + updated_path_planning_state[3];
    
    updated_path_planning_state[2] = fmax(next_d2s,0.0);
    
    return updated_path_planning_state;
}

vector<double> compute_next_ds(vector<double> path_planning_state){
    vector<double> updated_path_planning_state = compute_next_d2s(path_planning_state);

    double next_ds = path_planning_state[1] + updated_path_planning_state[2];
    std::cout<<next_ds<<endl;
    updated_path_planning_state[1] = next_ds;
    
    return updated_path_planning_state;
}

vector<double> compute_next_s(vector<double> path_planning_state){
    vector<double> updated_path_planning_state = compute_next_ds(path_planning_state);
    double next_s = path_planning_state[0] + updated_path_planning_state[1];
    
    updated_path_planning_state[0] = next_s;
    
    return updated_path_planning_state;
}


int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../../../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }
    
  vector<double> previous_path_ds;
  vector<double> previous_path_d2s;
  vector<double> previous_path_s;

    
  h.onMessage([&previous_path_s,&previous_path_ds,&previous_path_d2s,&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
      
      
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
            
        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];
            

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];
            
            
            // These variables are used for building the next set of waypoints
            double target_jerk = 0.00001;
            double max_acceleration = 0.004;
            double target_lane = 1;
            double target_velocity = 0.1;
            double car_acceleration = 0.0;
            double car_jerk = 0.0;
            
            int previous_size = previous_path_x.size();
            
            if(previous_size != 0)
            {
                car_acceleration = previous_path_d2s[50-previous_size];
                car_speed = previous_path_ds[50-previous_size];
                car_s = previous_path_s[50-previous_size];
            }
            
            previous_path_d2s.clear();
            previous_path_ds.clear();
            previous_path_s.clear();
            
            // The JSON message that will go back to the simulator
            json msgJson;
    
            vector<double> next_x_vals;
            vector<double> next_y_vals;
            
            vector<double> path_plannning_state;
            path_plannning_state.push_back(car_s);
            path_plannning_state.push_back(car_speed);
            path_plannning_state.push_back(car_acceleration);
            path_plannning_state.push_back(car_jerk);
            path_plannning_state.push_back(target_velocity);
            path_plannning_state.push_back(max_acceleration);
            path_plannning_state.push_back(target_jerk);
            
            
            vector<double> next_x_vals_smooth;
            vector<double> next_y_vals_smooth;
            
            for(int i=0;i<100;i++){
                
                path_plannning_state = compute_next_s(path_plannning_state);
                vector<double> next_xy = getXY(path_plannning_state[0], 6, map_waypoints_s, map_waypoints_x, map_waypoints_y);
                
                next_x_vals.push_back(next_xy[0]);
                next_y_vals.push_back(next_xy[1]);
                
                
                previous_path_s.push_back(path_plannning_state[0]);
                previous_path_ds.push_back(path_plannning_state[1]);
                previous_path_d2s.push_back(path_plannning_state[2]);
                
            }
            
            vector<double> ptsx;
            vector<double> ptsy;
            
            for(int i=0;i<next_x_vals.size();i+=25)
            {
                ptsx.push_back(next_x_vals[i]);
                ptsy.push_back(next_y_vals[i]);
            }
            
            // create a spline
            tk::spline s;
            s.set_points(ptsx, ptsy);
            
            for(int i=0;i<100;i++){
                double x_point = next_x_vals[i];
                double y_point = s(x_point);
                
                next_x_vals_smooth.push_back(x_point);
                next_y_vals_smooth.push_back(y_point);
                
            }
            
            msgJson["next_x"] = next_x_vals_smooth;
            msgJson["next_y"] = next_y_vals_smooth;
            
          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
          
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
