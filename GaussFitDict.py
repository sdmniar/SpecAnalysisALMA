# Rest Frequencies (in GHz)
freqs = {
    "CO1-0": 115.2712018,
    "CO2-1": 230.538,
    "CO3-2": 345.7959899,
    "CO4-3": 461.0407682,
    "CO5-4": 576.2679305,
    "CO6-5": 691.4730763,
    "CO7-6": 806.651806,
    "CO8-7": 921.7997,
    "CO9-8": 1036.912393,
    "CO10-9": 1151.985452,
    "CO11-10": 1267.014486,
    "CO12-11": 1381.995105,
    "CO13-12": 1496.922909,
    "CO14-13": 1611.793518,
    "CI1-0": 492.160651,
    "CI2-1": 809.34197,
    "CII": 1900.539,
}

# FITS files with PB correction
fits_list = {
    "BR1202": {
        "CI1-0":[
            "Data/BR/CI1-0/BR1202-0725_CI1-0_cleansub_spw1_o2.pbcor.image.fits",
        ],

        "CO5-4": [
            "Data/BR/CO5-4/BR1202-0725_CO5-4_cleansub_spw0_o1.pbcor.image.fits",
        ],

        "CO10-9": [
            "Data/BR/CO10-9/BRI1202-0725_CO10-9_cleansub_spw1_o1.pbcor.image.fits",
        ],

        "CO13-12":[
            "Data/BR/CO13-12/BR1202-0725_CO13-12_cleansub_spwALL_o1.pbcor.image.fits",
        ],

        "CO14-13": [
            "Data/BR/CO14-13/BR1202-0725_CO14-13_cleansub_spw0_o2.pbcor.image.fits",
        ],

        "CII": [
            "Data/BR/CII/BR1202-0725_CII_cleansub_spw04_o2.pbcor.image.fits"
        ],
    },

    "SGP38326": {

        "CO4-3": [
            "Data/SGP/CO4-3/SGP38326_CO4-3_cleansub_spw0_o1.pbcor.image.fits",
        ],

        "CI1-0": [
            "Data/SGP/CI1-0/SGP-38326_CI1-0_cleansub_spw5_o1.pbcor.image.fits",
        ],

        "CO5-4":[
            "Data/SGP/CO5-4/SGP38326_CO5-4_cleansub_spw3_o1.pbcor.image.fits",
        ],

        "CO7-6 & CI2-1": [
            "Data/SGP/CO7-6_CI2-1/SGP38326_CO7-6_CI2-1_cleansub_spw0_o1.pbcor.image.fits",
        ], 

        "CII": [
            "Data/SGP/CII/SGP38326_CII_cleansub_spwALL_o1.pbcor.image.fits",
        ],

    },
}

# .pb files
pb_list = {
    "BR1202": {
        "CI1-0":[
            "Data/BR/CI1-0/BR1202-0725_CI1-0_cleansub_spw1_o2.pb.fits",
        ],

        "CO5-4": [
            "Data/BR/CO5-4/BR1202-0725_CO5-4_cleansub_spw0_o1.pb.fits",
        ],

        "CO10-9": [
            "Data/BR/CO10-9/BRI1202-0725_CO10-9_cleansub_spw1_o1.pb.fits",
        ],

        "CO13-12":[
            "Data/BR/CO13-12/BR1202-0725_CO13-12_cleansub_spwALL_o1.pb.fits",
        ],

        "CO14-13": [
            "Data/BR/CO14-13/BR1202-0725_CO14-13_cleansub_spw0_o2.pb.fits",
        ],

        "CII": [
            "Data/BR/CII/BR1202-0725_CII_cleansub_spw04_o2.pb.fits"
        ],
    },

    "SGP38326": {

        "CO4-3": [
            "Data/SGP/CO4-3/SGP38326_CO4-3_cleansub_spw0_o1.pb.fits", 
        ],

        "CI1-0": [
            "Data/SGP/CI1-0/SGP-38326_CI1-0_cleansub_spw5_o1.pb.fits",
        ],

        "CO5-4":[
            "Data/SGP/CO5-4/SGP38326_CO5-4_cleansub_spw3_o1.pb.fits", 
        ],

        "CO7-6 & CI2-1": [
            "Data/SGP/CO7-6_CI2-1/SGP38326_CO7-6_CI2-1_cleansub_spw0_o1.pb.fits",
        ], 

        "CII": [
            "Data/SGP/CII/SGP38326_CII_cleansub_spwALL_o1.pb.fits",
        ],

    },
}

# FITS files without PB correction
noPBcor_list = {
    "BR1202": {
        "CI1-0":[
            "Data/BR/CI1-0/BR1202-0725_CI1-0_cleansub_spw1_o2.image.fits",
        ],

        "CO5-4": [
            "Data/BR/CO5-4/BR1202-0725_CO5-4_cleansub_spw0_o1.image.fits",
        ],

        "CO10-9": [
            "Data/BR/CO10-9/BRI1202-0725_CO10-9_cleansub_spw1_o1.image.fits",
        ],

        "CO13-12":[
            "Data/BR/CO13-12/BR1202-0725_CO13-12_cleansub_spwALL_o1.image.fits",
        ],

        "CO14-13": [
            "Data/BR/CO14-13/BR1202-0725_CO14-13_cleansub_spw0_o2.image.fits",
        ],

        "CII": [
            "Data/BR/CII/BR1202-0725_CII_cleansub_spw04_o2.image.fits"
        ],
    },

    "SGP38326": {

        "CO4-3": [
            "Data/SGP/CO4-3/SGP38326_CO4-3_cleansub_spw0_o1.image.fits", 
        ],

        "CI1-0": [
            "Data/SGP/CI1-0/SGP-38326_CI1-0_cleansub_spw5_o1.image.fits",
        ],

        "CO5-4":[
            "Data/SGP/CO5-4/SGP38326_CO5-4_cleansub_spw3_o1.image.fits",
        ],

        "CO7-6 & CI2-1": [
            "Data/SGP/CO7-6_CI2-1/SGP38326_CO7-6_CI2-1_cleansub_spw0_o1.image.fits",
        ], 

        "CII": [
            "Data/SGP/CII/SGP38326_CII_cleansub_spwALL_o1.image.fits",
        ],

    },
}