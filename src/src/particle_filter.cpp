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
#include <numeric>
#include "helper_functions.h"

using std::string;
using std::vector;
using std::normal_distribution;
void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  std::default_random_engine gen;
  num_particles = 100;  // TODO: Set the number of particles

  
  double std_x, std_y, std_theta;
  std_x = std[0];
  std_y = std[1];
  std_theta = std[2];
  normal_distribution<double> ini_dist_x(0,std_x);
  normal_distribution<double> ini_dist_y(0,std_y);
  normal_distribution<double> ini_dist_theta(0,std_theta);
  for (int i = 0;i<num_particles;++i){
    
    particles[i].id = i;
    particles[i].x = x;
    particles[i].y = y;
    particles[i].theta = theta;
    particles[i].weight = 1.0;
    
    particles[i].x += ini_dist_x(gen);
    particles[i].y += ini_dist_y(gen);
    particles[i].theta += ini_dist_theta(gen);
    
    
    
    
    //std::out<<"Sample"<<i+1<<" "<<sample_x<<" "<<sample_y<<" "<<sample_theta<<std::endl;
    
  }
  is_initialized = true;
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
  normal_distribution<double> dist_x(0,std_pos[0]);
  normal_distribution<double> dist_y(0,std_pos[1]);
  normal_distribution<double> dist_theta(0,std_pos[2]);
  for (unsigned int i=0;i<particles.size();++i){
    std::default_random_engine gen;
    if (fabs(yaw_rate)>0.00001){
      particles[i].x = particles[i].x + velocity*(sin(particles[i].theta + yaw_rate*delta_t)-sin(particles[i].theta))/yaw_rate;
      particles[i].y = particles[i].y + velocity*(cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t))/yaw_rate;
      particles[i].theta = particles[i].theta + yaw_rate*delta_t;
    } else {
      particles[i].x += velocity* delta_t * cos(particles[i].theta);
      particles[i].y += velocity* delta_t * sin(particles[i].theta);
    }
   
    
    particles[i].x += dist_x(gen);
    particles[i].y += dist_y(gen);
    particles[i].theta += dist_theta(gen);
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
  for (unsigned int i=0;i<observations.size();i++){
    LandmarkObs Obs = observations[i];
    double min_dist = std::numeric_limits<double>::max();
    int map_id = -1;
    for (unsigned int j=0;j<predicted.size();j++){
      LandmarkObs pre = predicted[j];
      double cur_dist = dist(Obs.x,Obs.y,pre.x,pre.y);
      if (cur_dist<min_dist){
        min_dist = cur_dist;
        map_id = pre.id;
      }
      observations[i].id = map_id;
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
  for(int i=0; i<num_particles;++i){
    
    double p_x = particles[i].x;
    double p_y = particles[i].y;
    double p_theta = particles[i].theta;
    
    vector<LandmarkObs> predictions;
    for (unsigned int j=0; j<map_landmarks.landmark_list.size();j++){
      double landmark_x = map_landmarks.landmark_list[j].x_f;
      double landmark_y = map_landmarks.landmark_list[j].y_f;
      int landmark_id = map_landmarks.landmark_list[j].id_i;
      
      if (pow(fabs(landmark_x-p_x),2) + pow(fabs(landmark_y-p_y),2) <= pow(sensor_range,2)){
        predictions.push_back(LandmarkObs{landmark_id, landmark_x, landmark_y});
      }
    }
    
    vector<LandmarkObs> transformed_obj;
    for (unsigned int k = 0; k<observations.size();k++){
      double transformed_x = cos(p_theta)*observations[k].x - sin(p_theta)*observations[k].y + p_x;
      double transformed_y = sin(p_theta)*observations[k].x + cos(p_theta)*observations[k].y + p_y;
      transformed_obj.push_back(LandmarkObs{observations[k].id, transformed_x, transformed_y});
    }
    
    dataAssociation(predictions, transformed_obj);


    particles[i].weight = 1.0;
    
    for (unsigned int m = 0; m<transformed_obj.size();m++){
      double mu_x, mu_y, pre_x, pre_y;
      mu_x = transformed_obj[m].x;
      mu_y = transformed_obj[m].y;
      
      int association_prediction = transformed_obj[m].id;
      
      for (unsigned int n = 0; n<predictions.size();n++){
        if (predictions[n].id == association_prediction){
          pre_x = predictions[n].x;
          pre_y = predictions[n].y;
        }
      }
      double sig_x = std_landmark[0];
      double sig_y = std_landmark[1];
      
      double up_w = (1/(2*M_PI*sig_x*sig_y)) * exp(-(pow(pre_x - mu_x,2)/(2*pow(sig_x, 2)) + (pow(pre_y-mu_y,2)/(2*pow(sig_y,2)))));
      
      particles[i].weight *= up_w;
    }
//    particles[i].weight = up_w;
//    weights[i] = particles[i].weight;
  }
  

}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  vector<Particle> new_particles(num_particles);
  std::default_random_engine gen;
  //vector<double> weights;
  for (int i = 0; i<num_particles;i++){
    weights.push_back(particles[i].weight);
  }
  
  std::uniform_int_distribution<int> uniintdist(0,num_particles-1);
  auto index = uniintdist(gen);
  
  double max_weight = *max_element(weights.begin(), weights.end());
  
  std::uniform_real_distribution<double> unirealdist(0.0, max_weight);
  
  double beta = 0.0;
  
  for (int i = 0; i < num_particles; i++){
    beta += unirealdist(gen)*2.0;
    while (beta > weights[index]){
      beta -= weights[index];
      index = (index+1)%num_particles;
    }
    new_particles.push_back(particles[index]);
  }
  particles = new_particles;
  //weights = particles.weights;
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