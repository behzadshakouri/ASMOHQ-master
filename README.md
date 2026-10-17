# ASMOHQ

## About

ASMOHQ is a C++/Qt research utility for generating activated-sludge process
models as OpenHydroQual projects. It includes model creators for fixed and
flexible configurations, example `.ohq` files, influent data, and result-grid
and VTK export helpers.

The generated models and their process parameters should be reviewed in
OpenHydroQual before they are used for scientific interpretation or design.

## OpenHydroQual resources

| Resource | Link |
|---|---|
| Project website | [openhydroqual.com](https://openhydroqual.com/) |
| Documentation and tutorials | [OpenHydroQual knowledge base](https://openhydroqual.com/knowledge-base/) |
| Searchable examples | [examples.openhydroqual.com](https://examples.openhydroqual.com/) |
| Example source repository | [ArashMassoudieh/OpenHydroQual-Examples](https://github.com/ArashMassoudieh/OpenHydroQual-Examples) |
| Source code and releases | [ArashMassoudieh/OpenHydroQual](https://github.com/ArashMassoudieh/OpenHydroQual) |
| Video tutorials | [OpenHydroQual on YouTube](https://www.youtube.com/@openhydroqual218) |

## Repository contents

- `ASMOHQ/`: application source and model-building components;
- `simple_model_files/`: compact model inputs and examples;
- `modelcreator*.cpp/.h`: fixed and flexible OpenHydroQual model generators;
- `CreatedModel*.ohq`: generated example projects;
- `DeNit_Influent_Lump.txt`: example influent data; and
- `resultgrid.*` and `vtk.h`: result handling and visualization support.

## Building

Open `ASMOHQ.pro` in a compatible Qt development environment, configure the
required local OpenHydroQual resources, and build the application. This is
research software; dependency versions and file paths may require adjustment
for the local system.

## Citation

If you use this repository, please cite the associated ASM1 surrogate-modeling
publication:

Shakouri, B., Murthy, S., De Clippeleir, H., Lesnik, K., & Massoudieh, A.
(2026). “Surrogate modeling of ASM1 using feedforward and LSTM neural networks.”
*Journal of Water Process Engineering, 81*, 109246.
[https://doi.org/10.1016/j.jwpe.2025.109246](https://doi.org/10.1016/j.jwpe.2025.109246)
· [Google Scholar record](https://scholar.google.com/citations?view_op=view_citation&hl=en&user=0lTtuMoAAAAJ&citation_for_view=0lTtuMoAAAAJ:4JMBOYKVnBMC)
