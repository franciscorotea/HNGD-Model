"""A module with helper functions used to write the output file from a 
simulation."""

def open_results_file(filename, node, unit_time='s', unit_temperature='K'):
    
    output_file = open(filename, 'w')
    
    output_file.write(f'Results for node = {node}\n')
    
    # Write header.
    
    match unit_time:
        case 's':
            time_str = 'time [s]'
        case 'm':
            time_str = 'time [m]'
        case 'h':
            time_str = 'time [h]'
        case 'd':
            time_str = 'time [d]'
            
    match unit_temperature:
        case 'K':
            temperature_str = 'temperature [K]'
        case 'C':
            temperature_str = 'temperature [°C]'
    
    headers = [time_str, temperature_str, 'growth [/]', 'nucleation [/]', 'Css [wt.ppm]', 'Cp [wt.ppm]']
    
    for header in headers:
        output_file.write(f'{header:^16}')
    
    output_file.write('\n')
    
    return output_file

def write_results(output_file, time, temperature, growth_fraction, nucleation_fraction, 
                  Css, Cp, unit_time='s', unit_temperature='K'):
    
    match unit_time:
        case 's':
            time = time
        case 'm':
            time = time/60
        case 'h':
            time = time/3600
        case 'd':
            time = time/86400
    
    match unit_temperature:
        case 'K':
            temperature = temperature
        case 'C':
            temperature = temperature - 273
    
    results = [time, temperature, growth_fraction, nucleation_fraction, Css, Cp]
    
    for result in results:
        output_file.write(f'{result:^16.4E}')
        
    output_file.write('\n')
    
def close_results_file(file):
    file.close()