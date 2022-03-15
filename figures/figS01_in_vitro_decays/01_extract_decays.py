from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import picoquant_tttr_sin_corr as pq

# This script (unlike most other scripts in the figure generation
# folder) directly accesses raw photon arrival time data (.SPT) in order
# to generate time-domain lifetime decays. Extracting from the raw data
# is kind of slow, so the output here is saved to csv and read in by the
# plotting python script in this same folder.

# generate some paths
current_dir = Path.cwd()
manuscript_path = current_dir.parents[1]
data_path = manuscript_path/'source_data'/'in_vitro_characterization'/'pH_temp_sample_ptu'
image_path_list = sorted(data_path.glob('*.ptu'))

# ndarray. pH, t (will drop y and x on parse)
pH_decays = np.zeros((8, 400))

# read in the images and extract decays
for j, file_path in enumerate(image_path_list):
    print("File %d, reading header..." % (j + 1), end='')
    tags = pq.parse_tttr_header(file_path, verbose=False)
    print("done")
    print("Num tags:", sum(len(tags[t]['values']) for t in tags))
    print("Num unique tags:", len(tags))
    print()
    # We loop over the frames one at a time:
    frames = pq.generate_picoharp_t3_frames(file_path, tags, verbose=False)
    for i, f in enumerate(frames):
        print("Parsing frame ", i, '... ', sep='', end='')
        parsed_frame = pq.parse_picoharp_t3_frame(
            records=f,
            tags=tags,
            verbose=True,
            show_plot=False,
            sin_corr_value=10, # equal to the scan speed / 10
            sinusoid_correction=True) 
        print("done.")
        trace_temp = pq.parsed_frame_to_histogram(parsed_frame,
                                                  x_pix_per_bin=512,
                                                  y_pix_per_bin=512,
                                                  t_pix_per_bin=4)
        np.add(pH_decays[j, :], trace_temp[:400, 1, :, :].reshape(400),
               out=pH_decays[j, :])

# save the output to txt
np.savetxt('2021-11-18_PBS_35C_decays.csv', pH_decays)


