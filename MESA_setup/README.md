# Wind Scheme Mapping

This JSON file **wind_schemes.json** helps you translate the wind_scheme ID numbers from your MESA history files into actual physical models and paper references.

When you run MESA with the run_star_extras.f90 provided in the Atlas, the code saves a number (like 24.0 or 45.5) in the wind_scheme column. This file acts as a dictionary so you can instantly see which paper that number refers to.

## What is inside

Each entry in the schemes dictionary uses the ID as a key and provides:

- Label: The short name of the recipe (e.g., GM23, V01)
- Author/year: Who wrote the paper and when
- Category: The regime of the wind (thin, thick, cool, LBV, etc.)
- full_citation: The complete bibliographic reference
- url: A link to the paper on ADS or arXiv

## How to use it

You can load this file in Python with just a few lines of code to label your plots or process your data:
  ```python
import json

with open('wind_scheme_map.json', 'r') as f:
    wind_map = json.load(f)

# (Note: MESA IDs are floats, so we convert to string to match the JSON keys)
# Example below
scheme_id = "24.0"
info = wind_map['schemes'].get(scheme_id)

if info:
    print(f"This is the {info['label']} model ({info['author']})")
    print(f"Category: {info['category']}")
  ```
This makes it much easier to keep track of which mass-loss recipe is active at different stages of the stellar evolution without having to check the Fortran source code every time.
