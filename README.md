# Stellar_Winds_Atlas
<img width="84" height="84" alt="image" src="https://github.com/user-attachments/assets/b61e70a4-7f3d-465c-99b5-c08470e20210" /> 

Code for the Stellar winds atlas paper for black hole formation. It uses the [StarEstate population synthesis tool](https://github.com/AmedeoRom/StarEstate) to generate entire populations of single stars for the Milky Way or elliptical galaxies.

Most of the data not included because too heavy to be uploaded. The only data included is the observation csv files and the final stellar population from one model
Use the GalPop.ipynb to produce a galactic population of massive stars. Use the following procedure:

1. In Section 1 create the samplers that will be used to quickly calculate the initial distribution of stars once the population synthesis will start (MODE = 'create_and_save_samplers'). It should take a few hours
2. After samplers are saved, the can be loaded without being calculated again by setting MODE = 'load_samplers_and_generate' with the same box you used to calculate samplers
3. By setting NUM_STARS, draw a population of how many stars you want (WATCH OUT: for a population that is too big, you might run out of RAM fast!) -- The Galactic position distribution is useless in this research. The code just calculates it out of consistency with the main StarEstate code
4. In Section 1.1, you can load your generated distribution of stars in case you already ran the code
5. In Section 1.2, use the code to match your stellar distributions with your MESA data files in terms of metallicity and ZAMS mass. Histograms on pre- and post-quantization should automatically appear. **Remember to separate your MESA folders firstly as a function of metallicity in solar masses in the format e.g. 0.1Zsun, 1Zsun, 2Zsun, etc., then in each metallicity folder put your masses in the format e.g. 0.1, 10, 200, etc.**
6. In Section 1.2.1 create the population from MESA simulations (what_to_do = 'generate'). You will create one or more parquet files (efficient at handling big datasets) that can be later loaded by setting what_to_do = 'just_extract' **REMEMBER TO SET THE RIGHT TRANSITION TO WR** 
7. In Section 1.3 define variables and functions, and load the observation samples
8. Plot the HRD. It is the same as in Figure 4 of the Stellar Winds Atlas I paper, but only with two panels. Here below what you would plot by using the parquet files that were provided

<img width="1692" height="1050" alt="image" src="https://github.com/user-attachments/assets/695179ca-950f-47af-a18c-ce5e69bcc76d" />
