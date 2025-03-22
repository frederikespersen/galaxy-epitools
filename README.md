# Galaxy EpiTools

This repository contains a set of custom tools for Galaxy Server. These tools are originally developed for EpiBiologics.

## Contents
The following tools, including both the Galaxy tool XML and function python scripts, are found in ``galaxy-epitools/``:
* ``ab-ngs-abundance``: A tool for counting distinct Ab sequences by their framework and CDRs, correcting for bias via metadata.
* ``calc-tada``: A tool for calculating [TA-DA Ab developability scores](https://www.tandfonline.com/doi/full/10.1080/19420862.2022.2080628#d1e434) from Schrodinger protein descriptors (See ``schrodinger-prot-desc``).
* ``calc-tap``: A tool for calculating [TAP Ab developability scores](https://doi.org/10.1073/pnas.1810576116) from ABodyBuilder2 models (See ``immunebuilder``).
* ``demultiplex-dorado``: A tool for demultiplexing Oxford Nanopore reads by the barcodes called by Dorado during basecalling.
* ``extract-ab-region``: A tool for extracting specific antibody regions according to their numbering by [ANARCI](https://github.com/oxpig/ANARCI).
* ``fasta-to-tabular``: A tool for extracting the header metadata in a FASTA file into a tabular format.
* ``immunebuilder``: A tool for running structure prediction with ABodyBuilder2 or NanoBodyBuilder2 from [ImmuneBuilder](https://github.com/oxpig/ImmuneBuilder).
* ``lima``: A small implementation of PacBio's [Lima](https://lima.how) for demultiplexing reads. This tool is developed with a very custom set of settings and is not intended for general use.
* ``paste``: An implementation of the basic shell command ``paste``.
* ``schrodinger-prot-desc``: A tool for calculating Ab protein descriptors with [Schrodinger](https://support.schrodinger.com/s/article/827119).
* ``sequence-distance``: A tool for determining a sequence distance matrix between two sets of sequences.
* ``sequence-liabilities``: A tool for finding sequence liabilities with regular expressions.
* ``space2``: A tool for epitope clustering antibodies by their structure with [SPACE2](https://github.com/oxpig/SPACE2).
* ``translate-sam``: A tool for translating reads mapped to a reference template, in the reading frame of the template.


## Getting started
Add these tools to your local Galaxy server instance in ``galaxy/tools/epitools/``
and add the following section to your ``galaxy/config/tools_conf.xml``:

```xml
<toolbox monitor="true">
  <!-- [... other tools ...] -->
  <section id="epitools" name="EpiTools">
    <tool file="epitools/ab-ngs-abundance/ab-ngs-abundance.xml"/>
    <tool file="epitools/calc-tada/calc-tada.xml"/>
    <tool file="epitools/calc-tap/calc-tap.xml"/>
    <tool file="epitools/demultiplex-dorado/demultiplex-dorado.xml"/>
    <tool file="epitools/extract-ab-regions/extract-ab-regions.xml"/>
    <tool file="epitools/fasta-to-tabular/fasta-to-tabular.xml"/>
    <tool file="epitools/immunebuilder/immunebuilder.xml"/>
    <tool file="epitools/lima/lima.xml"/>
    <tool file="epitools/paste/paste.xml"/>
    <tool file="epitools/schrodinger-prot-desc/schrodinger-prot-desc.xml"/>
    <tool file="epitools/sequence-distance/sequence-distance.xml"/>
    <tool file="epitools/sequence-liabilities/sequence-liabilities.xml"/>
    <tool file="epitools/space2/space2.xml"/>
    <tool file="epitools/translate-sam/translate-sam.xml"/>
  </section>
  <!-- [... other tools ...] -->
</toolbox>
```

Make sure to install any dependencies through the Galaxy server Admin interface.
