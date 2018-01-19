After you index and quantify you need to move the files into the following
folder structure so that it works with the R script. Or just fix the bloody R
script and kalipso scripts to output in a sensible order. Also see the
reference file used in the index for an example of the cDNA and ncRNA file, but
look this up to check (ie which to use cDNA, CDS, etc).

I also had issues with the biomart, if you get an error that your database is
missing just try to re-install the thing (see scrpt).
.
├── examples
├── index1.sh
├── my_meta.tsv
├── myResults
│   ├── CRI1
│   │   └── kalipso
│   │       ├── abundance.h5
│   │       ├── abundance.tsv
│   │       └── run_info.json
│   ├── CRI2
│   │   └── kalipso
│   │       ├── abundance.h5
│   │       ├── abundance.tsv
│   │       └── run_info.json
│   ├── CRI3
│   │   └── kalipso
│   │       ├── abundance.h5
│   │       ├── abundance.tsv
│   │       └── run_info.json
│   ├── CRI4
│   │   └── kalipso
│   │       ├── abundance.h5
│   │       ├── abundance.tsv
│   │       └── run_info.json
│   ├── KO1
│   │   └── kalipso
│   │       ├── abundance.h5
│   │       ├── abundance.tsv
│   │       └── run_info.json
│   ├── KO2
│   │   └── kalipso
│   │       ├── abundance.h5
│   │       ├── abundance.tsv
│   │       └── run_info.json
│   ├── KO3
│   │   └── kalipso
│   │       ├── abundance.h5
│   │       ├── abundance.tsv
│   │       └── run_info.json
│   └── KO4
│       └── kalipso
│           ├── abundance.h5
│           ├── abundance.tsv
│           └── run_info.json
├── quant2.sh
├── raw_reads
│   └── raw_reads -> /media/dwheeler/Seagate Expansion Drive/201712JHOL/201712JHOL/raw_reads/
├── README.md
├── ref
│   ├── Homo_sapiens.GRCh38.allrna.fa.gz
│   └── transcripts.idx
├── results
│   ├── SRR493366
│   │   └── kallisto
│   │       └── abundance.h5
│   ├── SRR493367
│   │   └── kallisto
│   │       └── abundance.h5
│   ├── SRR493368
│   │   └── kallisto
│   │       └── abundance.h5
│   ├── SRR493369
│   │   └── kallisto
│   │       └── abundance.h5
│   ├── SRR493370
│   │   └── kallisto
│   │       └── abundance.h5
│   └── SRR493371
│       └── kallisto
│           └── abundance.h5
├── slueth2.R
└── slueth.R

33 directories, 39 files
