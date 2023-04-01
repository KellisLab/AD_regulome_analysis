#!/usr/bin/env python3

import scanpy as sc
import numpy as np
import anndata
import pandas as pd
import os
import argparse

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("-b", "--bed", required=True)
