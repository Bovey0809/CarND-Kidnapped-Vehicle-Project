/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::normal_distribution;
using std::string;
using std::vector;
using namespace std;


void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  // cout<<"Start Init x:\n"<<x<<endl;
  num_particles = 100;  // TODO: Set the number of particles
  // num_particles is the private value for ParticleFilter.
  // The double std[] is the sigma pos for x, y, and theta.
  // This line creates a normal (Gaussian) distribution for x
  random_device rd;
  default_random_engine gen(rd());

  normal_distribution<double> dist_x(x, std[0]);
  // Create normal distributions for y and theta
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);
  // Initialize the vector size.
  particles.resize(num_particles);
  weights.resize(num_particles);
  // The weights need to fill in weight value later.
  for (int i = 0; i < num_particles; ++i)
  {
    particles[i].id = i;
    particles[i].x = dist_x(gen);
    particles[i].y = dist_y(gen);
    particles[i].theta = dist_theta(gen);
    particles[i].weight = 1.0;
  }
  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  
  default_random_engine gen;

  normal_distribution<double> dist_x(0, std_pos[0]);
  // Create normal distributions for y and theta
  normal_distribution<double> dist_y(0, std_pos[1]);
  normal_distribution<double> dist_theta(0, std_pos[2]);
  for(unsigned int i = 0; i < num_particles; i++)
  {
    double theta = particles[i].theta;
    if(abs(yaw_rate) < 0.00001){
      // Straight Line
      particles[i].x += velocity * cos(theta) * delta_t;
      particles[i].y += velocity * sin(theta) * delta_t;
    }
    else
    {
      particles[i].x += velocity / yaw_rate * (sin(theta + yaw_rate * delta_t) - sin(theta));
      particles[i].y += velocity / yaw_rate * (cos(theta) - cos(theta + yaw_rate * delta_t));
      particles[i].theta += yaw_rate * delta_t;
    }
    particles[i].x += dist_x(gen);
    particles[i].y += dist_y(gen);
    particles[i].theta += dist_theta(gen);
    // cout << "After Prediction: \n";
    // cout << "x: " << particles[i].x << " y: "<< particles[i].y << " theta: " << particles[i].theta << endl;
    // cout << "weight: " << particles[i].weight<<endl;
  }
  
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
  // Use the nearest neighbor algorithm to find out the nearest one predicted.
  // The landmarksObs struct contains id, x, and y.
  

  for(unsigned int i = 0; i < observations.size(); i++)
  {
    LandmarkObs obs = observations[i];
    double min_dist = 999999;
    double curr_dist;
    int closest_landmark = 0;
    for(unsigned int j = 0; j < predicted.size(); j++)
    {
      LandmarkObs pos = predicted[j];
      curr_dist = dist(pos.x, pos.y, obs.x, obs.y);
      if (curr_dist < min_dist)
      {
        min_dist = curr_dist;
        closest_landmark = j;
        
      }
    }
    observations[i].id = predicted[closest_landmark].id;
  }
  
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  for(unsigned int i = 0; i < num_particles; i++)
  {
    // cout<<"Before Update Weights"<<endl;
    double p_theta = particles[i].theta;
    double p_x = particles[i].x;
    double p_y = particles[i].y;
    
    // Use sensor range to get predictions from map_landmarks.
    vector<LandmarkObs> predictions;
    for (unsigned int j = 0; j < map_landmarks.landmark_list.size(); j++)
    {
      // sensor range = 50 m
      // The occupied area should be a circle. The origin is the particle itself.
      // We should find out all the map landmarks in the circle.
      // for each particles we only need the landmarks in the sensor range.
      Map::single_landmark_s landmark = map_landmarks.landmark_list[j];
      double landmark_dist = dist(landmark.x_f, landmark.y_f, p_x, p_y);
      // cout<<"landmark_dist: "<<landmark_dist<<endl;
      if (landmark_dist < sensor_range)
      {
        LandmarkObs pred = {landmark.id_i, landmark.x_f, landmark.y_f};
        // cout<<"prediction's x:"<<pred.x<<endl;
        predictions.push_back(pred);
      }
    }
    // transform the observations into the map coordinate.
    vector<LandmarkObs> transformed_os;
    for(unsigned int j = 0; j < observations.size(); j++)
    {
      double x_obs = observations[j].x;
      double y_obs = observations[j].y;
      
      double x_map = p_x + cos(p_theta) * x_obs - (sin(p_theta) * y_obs);
      double y_map = p_y + sin(p_theta) * x_obs + (cos(p_theta) * y_obs);
      
      LandmarkObs new_obs = {observations[j].id, x_map, y_map};
      transformed_os.push_back(new_obs);
    }

    
    //std::cout<<"Observations before association: \n"<<observations[0].id;
    dataAssociation(predictions, transformed_os);
    // Reinit weight
    particles[i].weight = 1.0;
    // update weights
    for(unsigned int j = 0; j < transformed_os.size(); j++)
    {
      // calculate weights
      double sig_x = std_landmark[0];
      double sig_y = std_landmark[1];
      double x_obs = transformed_os[j].x;
      double y_obs = transformed_os[j].y;
      int id_obs = transformed_os[j].id;
      double mu_x;
      double mu_y;
      for(unsigned int k = 0; k < predictions.size(); k++)
      {
        int id_pre = predictions[k].id;
        if (id_obs == id_pre) {
          mu_x = predictions[k].x;
          mu_y = predictions[k].y;
        }
      }
      // cout<<"mu_x: \n"<<mu_x<<endl;
      double update_weight = multiv_prob(sig_x, sig_y, x_obs, y_obs, mu_x, mu_y);
      // cout<<"update weight"<<update_weight<<endl;
      particles[i].weight *= update_weight;
      // cout<<"after update weight: "<<particles[i].weight<<endl;
    }
    weights.push_back(particles[i].weight);
  }
}

void ParticleFilter::resample() {

  vector<Particle> new_particles (num_particles);
  // cout<<"Resample:\n";
  random_device rd;
  mt19937 generator(rd());
  discrete_distribution<int> dist_weight(weights.begin(), weights.end());
  for(unsigned int  i = 0; i < num_particles; i++)
  {      
    new_particles[i] = particles[dist_weight(generator)];
  }
  particles = new_particles;
  
  weights.clear();
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

