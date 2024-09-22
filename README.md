# Conditional Probability of Surface Rupture - CPSR
Source code for the CPSR described in *Mammarella et al. (2024)*  

This repository contains a MATLAB script that calculates the conditional probability of surface rupture (CPSR) based on specific input values. The input values must be provided in a separate text file (`INPUT.txt`) located in the same directory. 

## Prerequisites  
- MATLAB (R2023a or higher)  
- A text file (`INPUT.txt`) with the required input values formatted correctly (see below).

## Input File Format  
The input file (`INPUT.txt`) must contain a single line of comma-separated values. Ensure that the values are provided in the exact order listed below.  

### Order of values in the input file:  

1. **MSR - Magnitude Scaling relation**: (`0` = L14_Interplate; `1` = L14_Intraplate; `2` = T17)

2. **SoF - Style of Fault**: (`3` = normal; `4` = reverse; `5` = strike-slip)

3. **HDD - Hypocentral Depth Ratio**: must be chosen consistent with the fault style:
   | Style of Fault (SoF) |    *Italy*    |  *California* | *Great Basin* |    *Taiwan*   |    *Japan*    | *New Zealand* |  *Aggregates* |
   |----------------------|---------------|---------------|---------------|---------------|---------------|---------------|---------------|
   | *if* SoF = 3  → HDD  |     `ITA_N`    |              |      `GB_N`   |               |               |               |     `AGG_N`   |
   | *if* SoF = 4  → HDD  |     `ITA_R`    |              |               |     `TAI_R`   |     `JAP_R`   |               |     `AGG_R`   |
   | *if* SoF = 5  → HDD  |               |     `CA_S`    |               |               |     `JAP_S`   |      `NZ_S`   |     `AGG_S`   |  

4. **dip_mu**: Mean fault dip angle calculated using a normal (Gaussian) distribution. (e.g, 40° for reverse fault)

5. **dip_sigma**: Standard deviation of *dip_mu* (e.g., 2)

6. **t_d**: Scaling parameter used to determine the width of the error in the fault angle. It is used to compute *eps_dip*, which represents the tolerated error around the mean (*dip_mu*). A larger value of t_d results in a wider error interval, allowing for more variability in the fault angles. (e.g., 2.5)

7. **Zs_mu**: Mean seismogenic depth based on a normal (Gaussian) distribution. (e.g., 14 km)

8. **Zs_sigma**: Standard deviation of *Zs_mu* (e.g., 2) 

9. **t_z**: Scaling parameter used to adjust the width of the error margin in seismogenic depth calculations.  It calculates *eps_Z*, which represents the tolerated error around the mean (*Zs_mu*). A larger value of t_z results in a wider error interval, allowing for more variability in the seismogenic depths. (e.g., 1)   


### Example:  
`0,4,AGG_R,40,2,2.5,14,2,1`  


## Running the Script

1. Open MATLAB and navigate to the directory containing the script and `INPUT.txt`.
2. **Option 1**: Run the script from the MATLAB command window:

   ```matlab
   run('CPSR.m')
   ```
   **Option 2**: Alternatively, open the script in MATLAB and click the Run button (green play button) in the editor.

___________________________________________________________
### Note: 
- Abbreviation in point 1: L14 refers to MSR by *[Leonard (2014)](https://doi.org/10.1785/0120140087)* for Interplate or Intraplate (stable continental region); T17 refert to MSR by *[Thingbaijam et al. (2017)](https://doi.org/10.1785/0120170017)*.

- String values listed in point 3 (e.g., ITA_N, GB_N, AGG_R, etc.) are associated with the corresponding numerical values of HDD derived from the research. For detailed information on these values, please refer to **Table 1** in *Mammarella et al. (2024) Supplementary Materials.*
___________________________________________________________
## License  
This project is licensed under the BSD 3-Clause License. See the [LICENSE](https://github.com/MammarellaLisa/CPSR/blob/main/LICENSE) file for details.

## Reference  
Leonard, M. (2014) Self-consistent earthquake fault-scaling relations: Update and extension to stable continental strike-slip faults. Bulletin of the Seismological Society of America, 104 (6): 2953–2965. doi:10.1785/0120140087.  
Thingbaijam, K.K.S., Mai, P.M. and Goda, K. (2017) New empirical earthquake source-scaling laws. Bulletin of the Seismological Society of America, 107 (5): 2225–2246. doi:10.1785/0120170017.  




