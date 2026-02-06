<div align="center">
 <img width="100" height="100" alt="image" src="https://github.com/user-attachments/assets/b61e70a4-7f3d-465c-99b5-c08470e20210" /> 
  <h1>Stellar Winds Atlas</h1>
  <p>
    <strong>A collection of stellar wind prescriptions for MESA.</strong>
  </p>
</div>

---

This repository contains the default [MESA](https://docs.mesastar.org/en/latest/index.html) (version 24.08.1) evolutionary setup for the Stellar Winds Atlas paper series

---

<div align="center">
 <img width="300" height="300" alt="ARGUS" src="https://github.com/user-attachments/assets/1be8ee13-4cee-4c06-be57-68b6901d6dc2" />
  <h2>Introducing Argus: The Stellar Winds Builder</h2>
  <h3><em>The Architect of the Wind Atlas</em></h3>
</div>

> **The Myth:** Argus was the master shipwright who constructed the *Argo* under the guidance of Athena for Jason and his Argonauts, who met the Hesperides, i.e. the daughters of Atlas.  
> **The Tool:** Argus is an AI-assisted companion that constructs the MESA inlist_project files for the user based on the guidelines of the Atlas papers.

> [!CAUTION]
> LLMs represent a huge environmental cost during training and usage, in terms of greenhouse gas emissions, energy consumption, and water usage for their cooling. [Companies may greenwash their operational emissions, while not being fully transparent about total lifecyle impacts](https://www.societybyte.swiss/en/2025/07/15/environmental-impact-of-large-language-models-green-or-polluting/). LLMs are very high cost-inefficient tools that use a lot of power to "guess" the right answer. Your brain uses orders of magnitude less power to reach to the same answer. Try to study the inlists on your own first, please!

> [!CAUTION]
> Please, always scrutinize the physical validity of your models. LLMs can produce wrong results! Also, this protocol is strictly a navigational tool for selecting initial metallicity and stellar wind configurations, issuing warnings where appropriate. It does not govern the full spectrum of initial conditions required for your simulations. Argus facilitates model construction but is by no means a substitute for human-based modelling and analysis. Do not treat your stellar model as a blackbox, please.

### What is Argus?
Configuring stellar winds in MESA involves complex variable mapping in `run_star_extras.f90`. **Argus** is a configuration protocol (powered by a strictly structured JSON schema) that allows Large Language Models (like ChatGPT, Claude, or Gemini) to function as an expert MESA "shipwright."

Instead of manually editing the `&controls` namelist, you provide the **Argus.json** (".json" like **Jason and the Argonauts**💡) to an LLM, and it builds the `inlist_project` for you based on simple natural language requests.

### How to Sail with Argus

To use the **Argus Builder**, you do not need to install packages. You simply need to boot the builder in your preferred AI chat interface.

1.  **Download the Context:**
    * Download `Argus.json` from this repository.
    * Download `inlist_project`.

2.  **The Prompt:**
    Upload both files to your LLM (it has been so far tested only with Gemini 1.5 Pro) and paste the following prompt:

    > "I am uploading the `Argus.json` (configuration logic) and `inlist_project` (template). You are now **Argus**, the Stellar Winds Builder. Please read the JSON to understand the wind regimes, warnings, and citations. Let' start! "

3.  **Build Your Model:**
    Ask Argus to build a specific wind scheme or maybe ask questions.
    * *Example:* "Configure a star using the 'Dutch' wind scheme."
    * *Example:* "I need a model for a massive star using the Krticka 2024 thin winds and Vink 2011 thick winds. Warn me if the metallicity is too low."
    * *Example:* "What are my options for thin winds?"
    * *Example:* "Where does this winds formula come from \[links and DOIs are included\]?"
    * *Example:* "Was this recipe derived from observations or simulations? What is the parameter space of the data from which it was fit?"

4.  **Receive the Code:**
    Argus will output the exact changes needed for your `&controls` block, ensuring no syntax errors and preserving your existing settings. Now you can combine it with run_star_extras and simulate!

### Why use Argus?
* **Safety:** It follows strict constraints defined in `Argus.json`, warning you if parameters (like metallicity) leave the physical regime of validity (e.g., $Z < 0.1 Z_{\odot}$).
* **Citations:** It automatically appends the correct references for the wind prescriptions you choose.
* **Consistency:** It ensures the `x_character_ctrl` flags map correctly to the `run_star_extras.f90` logic.

<!--Code for the Stellar winds atlas paper for black hole formation. **It includes the adopted MESA setup and the setup for creating the OB/WR galactic population**

It uses the [StarEstate population synthesis tool](https://github.com/AmedeoRom/StarEstate) to generate entire populations of single stars for the Milky Way or elliptical galaxies (the galaxy shape is irrelevant in this paper).

Most of the data not included because too heavy to be uploaded. The only data included is the observation csv files and the final stellar population from one model
Use the GalPop.ipynb to produce a galactic population of massive stars. Use the following procedure:

1. In Section 1 create the samplers that will be used to quickly calculate the initial distribution of stars once the population synthesis will start (MODE = 'create_and_save_samplers'). It should take a few hours
2. After samplers are saved, the can be loaded without being calculated again by setting MODE = 'load_samplers_and_generate' with the same box you used to calculate samplers
3. By setting NUM_STARS, draw a population of how many stars you want (WATCH OUT: for a population that is too big, you might run out of RAM fast!) -- The Galactic position distribution is useless in this research. The code just calculates it out of consistency with the main StarEstate code
4. In Section 1.1, you can load your generated distribution of stars in case you already ran the code
5. In Section 1.2, use the code to match your stellar distributions with your MESA data files in terms of metallicity and ZAMS mass. Histograms on pre- and post-quantization should automatically appear. **Remember to separate your MESA folders firstly as a function of metallicity in solar masses in the format e.g. 0.1Zsun, 1Zsun, 2Zsun, etc., then in each metallicity folder put your masses in the format e.g. 0.1, 10, 200, etc.**
6. In Section 1.2.1 create the population from MESA simulations (what_to_do = 'generate'). You will create one or more parquet files (efficient at handling big datasets) that can be later loaded by setting what_to_do = 'just_extract'. The partition to multiple parquet files is used to not kill your RAM in case your population is huge! **REMEMBER TO SET THE RIGHT TRANSITION TO WR** 
7. In Section 1.3 define variables and functions, and load the observation samples
8. Plot the HRD. It is the same as in Figure 4 of the Stellar Winds Atlas I paper, but only with two panels. Here below what you would plot by using the parquet files that were provided

<img width="1692" height="1050" alt="image" src="https://github.com/user-attachments/assets/695179ca-950f-47af-a18c-ce5e69bcc76d" /> -->

---

### The Hall of Fame

A project like this doesn't sail alone. Special thanks to:

* **Marco Ruscone**: For teaching me how to use `.json` (the *Jason* to my Argonauts) to guide LLMs.
* **Lucas de Sá**: For suggesting the name **Argus** upon this shipwright.
* **My own subconscious**: For quite literally dreaming up the implementation logic into MESA one night, apparently deciding that astrophysics was more important than a good night's sleep.
