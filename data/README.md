# Data directory
Place your data (`.fcs`) in this directory. Make one directory for each sample/analyte and place all raw data files in the respective folder. For the downstream analysis to work you should have the same number of files (i.e. the same residue positions) for each analyte.

## Example
Say that you have two analytes, ab1 and ab2, for which you have tested 5 residue positions. You also have two positive and two negative controls. Your target is called Tar and you have tested the following positions:
- TAR312T
- TAR317Y
- TAR332G
- TAR343Y
- TAR354S
Make two directories, one named ab1 and one named ab2. In each directory you place all `.fcs` files beloning to that sample. In the end you should have the following file structure:
```
data/
├── AB1
│   ├── NEG_AB1(1).fcs
│   ├── NEG_AB1(2).fcs
│   ├── POS_AB1(1).fcs
│   ├── POS_AB1(2).fcs
│   ├── TAR312T_AB1.fcs
│   ├── TAR317Y_AB1.fcs
│   ├── TAR332G_AB1.fcs
│   ├── TAR343Y_AB1.fcs
│   └── TAR354S_AB1.fcs
├── AB1
│   ├── NEG_AB2(1).fcs
│   ├── NEG_AB2(2).fcs
│   ├── POS_AB2(1).fcs
│   ├── POS_AB2(2).fcs
│   ├── TAR312T_AB2.fcs
│   ├── TAR317Y_AB2.fcs
│   ├── TAR332G_AB2.fcs
│   ├── TAR343Y_AB2.fcs
│   └── TAR354S_AB2.fcs
└── README.md
```

Following this naming of the files will make everything look good. After running the analysis once for an analyte that subfolder will also contain a `defaults.toml` file which contains your run settings. Unless you use the `-d` tag when rerunning it will reuse thsoe settings. 