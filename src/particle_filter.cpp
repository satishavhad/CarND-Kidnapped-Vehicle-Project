/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>
#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).


	// To add Gaussian noise to each particle (x, y and theta).
	default_random_engine rand_gen;
	normal_distribution<double> dst_x(x, std[0]);
	normal_distribution<double> dst_y(y, std[1]);
	normal_distribution<double> dst_theta(theta, std[2]);
	
	// Number of particles to draw
	num_particles = 100;

	// Flag, if filter is initialized
	is_initialized = true;

	// Resize Vector of particles and weights of all particles to number of particles
	weights.resize(num_particles);
	particles.resize(num_particles);

	//Initialize particles
	for (int i = 0; i<num_particles; i++)
	{
		particles[i].id = i;
		particles[i].x = dst_x(rand_gen);
		particles[i].y = dst_y(rand_gen);
		particles[i].theta = dst_theta(rand_gen);
		particles[i].weight = 1.0;
	}

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	default_random_engine rand_gen;
	normal_distribution<double> dst_x(0, std_pos[0]);
	normal_distribution<double> dst_y(0, std_pos[1]);
	normal_distribution<double> dst_yaw(0, std_pos[2]);

	for (int i = 0; i<num_particles; i++)
	{
		//Handle division by zero case
		if (fabs(yaw_rate) > 0.001) {
			particles[i].x += velocity / yaw_rate * (sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta));
			particles[i].y += velocity / yaw_rate * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t));
		}
		else {
			particles[i].x += velocity*delta_t*cos(particles[i].theta);
			particles[i].y += velocity*delta_t*sin(particles[i].theta);
		}

		particles[i].theta += yaw_rate*delta_t;

		// Adding noise to x, y and theta using std::normal_distribution and std::default_random_engine
		particles[i].x += dst_x(rand_gen);
		particles[i].y += dst_y(rand_gen);
		particles[i].theta += dst_yaw(rand_gen);
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	
	double min_dst;
	int map_idfc;
	double current_dist, deltaX, deltaY;

	for (int j = 0; j < observations.size(); j++)
	{
		LandmarkObs obsrv = observations[j];

		//assign minimum distance to maximum value to begin with
		min_dst = numeric_limits<double>::max();
		
		//identification of the landmark id initialization
		map_idfc = -1;

		for (int i = 0; i < predicted.size(); ++i)
		{
			LandmarkObs landmark_pred = predicted[i];
			
			deltaX = (landmark_pred.x - obsrv.x);
			deltaY = (landmark_pred.y - obsrv.y);
			current_dist = deltaX*deltaX + deltaY*deltaY;
			
			if (current_dist < min_dst)
			{
				min_dst = current_dist;
				map_idfc = i;
			}
		}
		// set the observation ID to nearest neighbour
		observations[j].id = map_idfc; 
	}

}


void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

	for (int j = 0; j < num_particles; j++)
	{
		//Collect Landmark Predictions if map landmark locations predicted are within sensor range of the particle
		std::vector<LandmarkObs> landmark_predictions;
		for (auto landmark : map_landmarks.landmark_list)
		{
			LandmarkObs landmark_pred;
			landmark_pred.id = landmark.id_i;
			landmark_pred.x = landmark.x_f;
			landmark_pred.y = landmark.y_f;
			
			//Collect Landmark Predictions if map landmark locations predicted are within sensor range of the particle
			if (fabs(landmark_pred.x - particles[j].x)<= sensor_range && fabs(landmark_pred.y - particles[j].y) <= sensor_range)
				landmark_predictions.push_back(landmark_pred);
		}
		//Collected Landmark Predictions


		//observations to be transformed from vehicle to map coordinates
		std::vector<LandmarkObs> trnsf_obsrv;
		for (LandmarkObs obs_landmark : observations)
		{
			LandmarkObs tr_ob;
			tr_ob.id = obs_landmark.id;
			tr_ob.x = particles[j].x + obs_landmark.x * cos(particles[j].theta) - obs_landmark.y * sin(particles[j].theta);
			tr_ob.y = particles[j].y + obs_landmark.x * sin(particles[j].theta) + obs_landmark.y * cos(particles[j].theta);
			
			trnsf_obsrv.push_back(tr_ob);
		}
		//observations transformeded


		// call dataAssociation defined above
		dataAssociation(landmark_predictions, trnsf_obsrv);

		//Calculate Weights
		double total_probability = 1.0;
		for (int i = 0; i < trnsf_obsrv.size(); i++)
		{
			LandmarkObs observation_temp = trnsf_obsrv[i];
			LandmarkObs landmark_temp = landmark_predictions[observation_temp.id];

			double deltaX = (observation_temp.x - landmark_temp.x);
			double deltaY = (observation_temp.y - landmark_temp.y);
			double posterior_distr = exp(-(deltaX*deltaX / (2 * std_landmark[0] * std_landmark[0]) + deltaY*deltaY / (2 * std_landmark[1] * std_landmark[1]))) / (2.0*M_PI*std_landmark[0] * std_landmark[1]);
						
			total_probability *= posterior_distr;
		}

		particles[j].weight = total_probability;
		weights[j] = total_probability;
		//Weights are Calculated
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	std::discrete_distribution<int> index(weights.begin(), weights.end());
	std::vector<Particle> Resample_particles;

	default_random_engine rand_gen;

	for (int i = 0; i < num_particles; i++)
	{
		int temp = index(rand_gen);
		Resample_particles.push_back(particles[temp]);
	}
	particles = Resample_particles;
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
