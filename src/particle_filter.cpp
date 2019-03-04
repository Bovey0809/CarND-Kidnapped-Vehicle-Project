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

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  std::cout<<"Start Init \n";
  num_particles = 100;  // TODO: Set the number of particles
  // num_particles is the private value for ParticleFilter.
  // The double std[] is the sigma pos for x, y, and theta.
  // This line creates a normal (Gaussian) distribution for x
  float std_x = std[0];
  float std_y = std[1];
  float std_theta = std[2];
  std::default_random_engine gen;

  normal_distribution<double> dist_x(x, std_x);
  // Create normal distributions for y and theta
  normal_distribution<double> dist_y(y, std_y);
  normal_distribution<double> dist_theta(theta, std_theta);
  // Initialize the vector size.
  particles.resize(num_particles);
  weights.resize(num_particles);
  // The weights need to fill in weight value later.
  for (int i = 0; i < num_particles; ++i)
  {
    double sample_x, sample_y, sample_theta;

    //   where "gen" is the random engine initialized earlier.
    sample_x = dist_x(gen);
    sample_y = dist_y(gen);
    sample_theta = dist_theta(gen);
    particles[i].id = i;
    particles[i].x = sample_x;
    particles[i].y = sample_y;
    particles[i].theta = sample_theta;
    particles[i].weight = 1.0;
  }
  is_initialized = true;
  std::cout<<"New Particle size: "<<particles.size()<<"\n Weight: "<<particles[0].weight;
  std::cout<<"\n Init Done";
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  // Why should I add random Gaussian noise?
  float std_x = std_pos[0];
  float std_y = std_pos[1];
  float std_theta = std_pos[2];
  std::default_random_engine gen;

  normal_distribution<double> dist_x(0, std_x);
  // Create normal distributions for y and theta
  normal_distribution<double> dist_y(0, std_y);
  normal_distribution<double> dist_theta(0, std_theta);
  for(unsigned int i = 0; i < num_particles; i++)
  {
    // std::cout<<"Previous Particle: \n"<<particles[i].x<<"\n";
    // Particle *particle = &particles[i];
    double theta = particles[i].theta;
    if(fabs(yaw_rate) < 0.00001){
      // Straight Line
      particles[i].x += velocity * cos(theta) * delta_t + dist_x(gen);
      particles[i].y += velocity * sin(theta) * delta_t + dist_y(gen);
    }
    else
    {
      particles[i].theta += yaw_rate * delta_t;
      particles[i].x += velocity / yaw_rate * (sin(particles[i].theta) - sin(theta));
      particles[i].y += velocity / yaw_rate * (cos(theta) - cos(particles[i].theta));
    }    
    // std::cout << "Now Particle: \n"<< particles[i].x << "\n";
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
  double min_dist = 999999;

  for(unsigned int i = 0; i < observations.size(); i++)
  {
    LandmarkObs obs = observations[i];
    for(unsigned int j = 0; j < predicted.size(); j++)
    {
      LandmarkObs pos = predicted[j];
      double cur_dist = dist(pos.x, pos.y, obs.x, obs.y);
      if (cur_dist < min_dist)
      {
        observations[i].id = pos.id;
        min_dist = cur_dist;
      };
    }
  }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
  for(unsigned int i = 0; i < num_particles; i++)
  {
    double theta = particles[i].theta;

    // transform the observations into the map coordinate.
    vector<LandmarkObs> transformed_os;
    for(unsigned int j = 0; j < observations.size(); j++)
    {
      double x_obs = observations[i].x;
      double y_obs = observations[i].y;
      
      double x_map = particles[i].x + (cos(theta) * x_obs) - (sin(theta) * y_obs);
      double y_map = particles[i].y + (sin(theta) * x_obs) + (cos(theta) * y_obs);
      
      LandmarkObs new_obs = {observations[j].id, x_map, y_map};
      transformed_os.push_back(new_obs);
    }

    // Use sensor range to get predictions from map_landmarks.
    vector<LandmarkObs> predictions;
    for(unsigned int j = 0; j < map_landmarks.landmark_list.size(); j++)
    {
      // sensor range = 50 m
      // The occupied area should be a circle. The origin is the particle itself.
      // We should find out all the map landmarks in the circle.
      // for each particles we only need the landmarks in the sensor range.
      Map::single_landmark_s landmark = map_landmarks.landmark_list[j];
      double landmark_dist = dist(landmark.x_f, landmark.y_f, particles[i].x, particles[i].y);
      if (landmark_dist < sensor_range) { 
        LandmarkObs pred = {landmark.id_i, landmark.x_f, landmark.y_f};
        predictions.push_back(pred);
      }
    }
    //std::cout<<"Observations before association: \n"<<observations[0].id;
    dataAssociation(predictions, transformed_os);
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
      double update_weight = multiv_prob(sig_x, sig_y, x_obs, y_obs, mu_x, mu_y);
      particles[i].weight *= update_weight;
    }
  }
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
std::default_random_engine generator;

for(unsigned int  i = 0; i < num_particles; i++)
{
  std::discrete_distribution<int> dist(weights.begin(), weights.end());
  particles[i] = particles[dist(generator)];
}
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

