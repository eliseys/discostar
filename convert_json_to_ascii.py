import json



data = json.load(open('discostar/DATA/data.json'))


output_file = "HZ_Her-B.data"

with open(output_file, "w") as f:
    for n in data:
        JD, flux, orbit, prec = [], [], [], []

        JD = data[n]['B']['JD']
        flux = data[n]['B']['flux']
        orbit = data[n]['B']['orbit']
        prec = data[n]['B']['prec']

        for i, t in enumerate(JD):
            f.write("{}\t{}\t{}\t{}\n".format(t, flux[i], orbit[i], prec[i]))


output_file = "HZ_Her-V.data"

with open(output_file, "w") as f:
    for n in data:
        JD, flux, orbit, prec = [], [], [], []

        JD = data[n]['V']['JD']
        flux = data[n]['V']['flux']
        orbit = data[n]['V']['orbit']
        prec = data[n]['V']['prec']

        for i, t in enumerate(JD):
            f.write("{}\t{}\t{}\t{}\n".format(t, flux[i], orbit[i], prec[i]))

