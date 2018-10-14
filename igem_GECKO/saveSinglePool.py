from geckopy import GeckoModel
import cobra, pandas, os

model = GeckoModel('single-pool')

cobra.io.save_matlab_model(model, "ec_model_singlepool.mat")