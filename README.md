# radqm9_pipeline
This is a pipeline for formatting the radqm9 dataset for MACE. The pipeline can create different dataset configurations based on the properties within radqm9.

1. How is the data cleaned?
The data begins as a collection of json files. These files are read and merged into a single list. There are assumptions made on the format of a single data point:
- Each point is a python dictionary with the keys serving as fields for storing various data values. 
- The mol_id (one of the fields (key) within the data point dictionary) can be thought of as a "root". Data that have this mol_id correspond to the same molecule/optimization problem; however, each data point represents a unique charge/spin pair. 

Duplicates charge/spin pairs can happen, but are not considered "copies". Rather we need to check if it is an optimization that stopped and restarted. To do so, we check to see if an order can be established for the duplicate charge/spin pairs, meaning the end value of one should be the beginning of another (positions field). If no order can be defined, the duplicate pairs are all removed. 

After resolving duplicates, the data is downsampled as follows:
- For a charge of "0" and spin of "1", i.e a charge/spin pair of '0,1', we choose the first and the last points for geometries, energies, forces, resp, and mulliken. 
- for the remaining data we calculate the average force at each optimization snapshot. Using this we take the first point, the last point, and the point with the highest average force from the remaining optimization snapshots. 

Once downsampled, we need to search for data with "broken graphs" During the optimization, the molcule can split and form multiple non-connected objects. These are removed.

2. How is the data split?
The data, by default, is split 65/25/10 (train/val/test respectively).

We want to ensure equal representation of molecule types and prevent any form of data leakage. For example, not every data point has Flourine, but we want to make sure that Flourine is in the train/val/test as equally as possible. We also want to ensure molecule size (we define this by the sum of the atomic masses) is also equally represented. Even though we want to ensure some form of equal representation, we don't want mol_ids to be spread across splits. This means that all charge/spin pairs associated with the mol_id are fully contained in either train/val/test such.

To accomplish this we create a dictionary where a key is the sum of the atomic mass for each element in the molecule. Each key will be associated with a list of mol_ids that have that molecule sum. We can then take the defined percentage splits for each mass bin (key).

3. How is the data chunked?
Depending on the model, training on the entire dataset can be compute time intensive, requiring resources that are not available to every member of the community. Developers would also wish to explore ideas with rapid prototype datasets, but also investigate model performance based on dataset sizes. As a result, we present data chunks for the training set: 5%, 10%, 25%, 50%, 75%, and of course 100%.

Much like creating the train/val/test split, there are some desirable properties we wish to uphold within our chunks. 
- Each subsequent chunk contains the prior chunks, i.e., the 25% chunk contains the 10% and 5%, the 10% contains the 5%, etc. 
- We also want to ensure equal representation within these chunks. We follow a similar methods to the data split section above. 