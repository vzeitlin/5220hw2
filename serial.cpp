#include "common.h"
#include <cmath>
#include <vector>
#include <iostream>
#include <utility>

// Apply the force from neighbor to particle
void apply_force(particle_t& particle, particle_t& neighbor) {
    // Calculate Distance
    double dx = neighbor.x - particle.x;
    double dy = neighbor.y - particle.y;
    double r2 = dx * dx + dy * dy;

    // Check if the two particles should interact
    if (r2 > cutoff * cutoff)
        return;

    r2 = fmax(r2, min_r * min_r);
    double r = sqrt(r2);

    // Very simple short-range repulsive force
    double coef = (1 - cutoff / r) / r2 / mass;
    particle.ax += coef * dx;
    particle.ay += coef * dy;
}

// Integrate the ODE
void move(particle_t& p, double size) {
    // Slightly simplified Velocity Verlet integration
    // Conserves energy better than explicit Euler method
    p.vx += p.ax * dt;
    p.vy += p.ay * dt;
    p.x += p.vx * dt;
    p.y += p.vy * dt;

    // Bounce from walls
    while (p.x < 0 || p.x > size) {
        p.x = p.x < 0 ? -p.x : 2 * size - p.x;
        p.vx = -p.vx;
    }

    while (p.y < 0 || p.y > size) {
        p.y = p.y < 0 ? -p.y : 2 * size - p.y;
        p.vy = -p.vy;
    }
}

const float bin_size = 0.03;
static std::vector<std::pair<int, std::vector<int>>> bins;
static int num_bins;

void init_simulation(particle_t* parts, int num_parts, double size) {
	// You can use this space to initialize static, global data objects
    // that you may need. This function will be called once before the
    // algorithm begins. Do not do any particle simulation here
    int bin_div = (int) (size/bin_size);
    num_bins = (bin_div == size/bin_size) ? bin_div : bin_div + 1;

    //std::cout << "num bins: " << num_bins << " size: " << size << std::endl;

    bins = std::vector<std::pair<int,std::vector<int>>>(num_bins*num_bins, std::pair<int,std::vector<int>>(0, std::vector<int>()));

    
}

inline void compute_inter_bin(particle_t* parts, int bin_i1, int bin_i2){
    //std::cout << "(" << bin_i1 % num_bins << "," << bin_i1 / num_bins << ") <-> (" 
    //<< bin_i2 % num_bins << "," << bin_i2 / num_bins << ")\n";

    for(int i = 0; i < std::get<0>(bins[bin_i1]); i++){
        int curr_part =  std::get<1>(bins[bin_i1])[i];
        for(int j = 0; j < std::get<0>(bins[bin_i2]); j++){
            apply_force(parts[curr_part], parts[(std::get<1>(bins[bin_i2])[j])]);
        }
    }
}


void simulate_one_step(particle_t* parts, int num_parts, double size) {
    // Compute Forces
    /*
    for (int i = 0; i < num_parts; ++i) {
        parts[i].ax = parts[i].ay = 0;
        for (int j = 0; j < num_parts; ++j) {
            apply_force(parts[i], parts[j]);
        }
    }
    
    */
   for(int i = 0; i < num_parts; i++){
        //std::cout << "P: " << parts[i].x << "," << parts[i].y << " to bin " << (int)(parts[i].x / bin_size) << ","
        //<< (int)(parts[i].y / bin_size) << std::endl;
        parts[i].ax = parts[i].ay = 0;
        int curr_idx = (int)(parts[i].x / bin_size) + ((int)(parts[i].y / bin_size))*num_bins;
        if(++(std::get<0>(bins[curr_idx])) > std::get<1>(bins[curr_idx]).size()){
            std::get<1>(bins[curr_idx]).push_back(i);
        } else{
           std::get<1>(bins[curr_idx])[std::get<0>(bins[curr_idx])-1] = i;
        }

    }
    
    for(int i = 0; i < num_bins*num_bins; i++){
        bool x_bound = i % num_bins == num_bins - 1; 
        bool y_bound = i / num_bins == num_bins - 1;
        bool x_bound_top = i % num_bins == 0;
        bool y_bound_top = i / num_bins == 0;

        int y_next = i + num_bins;
        int y_prev = i - num_bins;
        
        if(!x_bound && !y_bound_top){
            compute_inter_bin(parts, i, y_prev + 1);
        }
        if(!x_bound && !y_bound){
            compute_inter_bin(parts, i, y_next + 1);
        }
        if(!x_bound){
            compute_inter_bin(parts, i, i + 1);
        }
        if(!x_bound_top){
            compute_inter_bin(parts, i, i - 1);
        }

        if(!x_bound_top && !y_bound_top){
            compute_inter_bin(parts, i, y_prev - 1);
        }
        if(!x_bound_top && !y_bound){
            compute_inter_bin(parts, i, y_next - 1);
        }
        if(!y_bound){
            compute_inter_bin(parts, i, y_next);
        }
        if(!y_bound_top){
            compute_inter_bin(parts, i, y_prev);
        }

        compute_inter_bin(parts, i, i);
    }
    

    for(int i = 0; i < num_bins*num_bins; i++){
        std::get<0>(bins[i]) = 0;
    }
    // Move Particles
    for (int i = 0; i < num_parts; ++i) {
        move(parts[i], size);
    }
    //std::cout << "\n\n";
}
