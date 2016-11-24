import pylab, os

x = []
y = []

for file in os.listdir(os.getcwd()):
    if file.endswith(".txt"):
        _file = open(file)
        for line in _file:
            x.append(float(line.split('\t')[0]))
            y.append(float(line.split('\t')[1]))
        _file.close()
        pylab.plot(x, y)
pylab.show()
