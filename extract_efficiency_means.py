#!/usr/bin/env python
"""
Extract mean efficiencies from validation histograms
"""
from __future__ import print_function

import ROOT
import sys

def extract_means(root_file_path):
    """Extract mean values from efficiency histograms"""
    
    # Open the ROOT file
    f = ROOT.TFile.Open(root_file_path)
    if not f or f.IsZombie():
        print("Error: Could not open file {0}".format(root_file_path))
        sys.exit(1)
    
    # Define the histogram names
    particles = {
        'Gluino': [1000, 1800, 2200],
        'Stop': [1000, 1800, 2200],
        'Stau': [550, 850, 1200]
    }
    
    selection_stages = ['tight']
    # selection_stages = ['trigPres', 'selection', 'tight']
    
    # Store results
    results = {}
    
    # Loop through all histograms
    for particle, masses in particles.items():
        results[particle] = {}
        for mass in masses:
            results[particle][mass] = {}
            for stage in selection_stages:
                # Construct histogram name
                hist_name = "h_{0}{1}_efficiency_{2}".format(particle.lower(), mass, stage)
                
                # Get histogram
                hist = f.Get(hist_name)
                if hist:
                    mean = hist.GetMean()
                    std = hist.GetStdDev()
                    entries = hist.GetEntries()
                    results[particle][mass][stage] = {
                        'mean': mean,
                        'std': std,
                        'entries': entries
                    }
                else:
                    print("Warning: Could not find histogram {0}".format(hist_name))
    
    f.Close()
    
    return results


def print_results(results):
    """Print results in a formatted table"""
    
    print("\n" + "="*80)
    print("EFFICIENCY MEANS SUMMARY")
    print("="*80)
    
    for particle, masses in results.items():
        print("\n{0}:".format(particle))
        print("-" * 80)
        print("{0:<12} {1:<20} {2:<12} {3:<12} {4:<10}".format('Mass (GeV)', 'Stage', 'Mean', 'Std Dev', 'Entries'))
        print("-" * 80)
        
        for mass, stages in sorted(masses.items()):
            for stage, values in stages.items():
                stage_name = {
                    # 'trigPres': 'Trigger+Prescale',
                    # 'selection': 'Trig+Pres+Selection',
                    'tight': 'Trig+Pres+Sel+Tight'
                }[stage]
                
                print("{0:<12} {1:<20} {2:<12.4f} {3:<12.4f} {4:<10.0f}".format(
                    mass, stage_name, values['mean'], values['std'], values['entries']))


def save_to_file(results, output_file):
    """Save results to a text file"""
    
    with open(output_file, 'w') as f:
        f.write("Particle,Mass_GeV,SelectionStage,Mean,StdDev,Entries\n")
        
        for particle, masses in results.items():
            for mass, stages in sorted(masses.items()):
                for stage, values in stages.items():
                    stage_name = {
                        # 'trigPres': 'TriggerPrescale',
                        # 'selection': 'TrigPresSelection',
                        'tight': 'TrigPresSelectionTight'
                    }[stage]
                    
                    f.write("{0},{1},{2},{3:.6f},{4:.6f},{5:.0f}\n".format(
                        particle, mass, stage_name, values['mean'], values['std'], values['entries']))


if __name__ == "__main__":
    # Check command line arguments
    if len(sys.argv) < 2:
        print("Usage: python extract_efficiency_means.py <root_file> [output.csv]")
        sys.exit(1)
    
    root_file = sys.argv[1]
    
    # Extract means
    results = extract_means(root_file)
    
    # Print results
    print_results(results)
    
    # Save to file if requested
    if len(sys.argv) >= 3:
        output_file = sys.argv[2]
        save_to_file(results, output_file)
        print("\nResults saved to {0}".format(output_file))
