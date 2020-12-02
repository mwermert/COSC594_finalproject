import csv
import sys
import matplotlib.pyplot as plt

markers = ['.', 'o', 'v', '^', '<', '>', '1', '2', '3', '4']

def main(argv):
    with open(argv[0]) as f:

        organismList = {}

        reader = csv.reader(f)
        count = 0
        for row in reader:
            if count > 0 and row[2] != '':
                # do something here with `row`
                if not row[2] in organismList:
                    organismList[row[2]] = []
                organismList[row[2]].append((row[1], row[3]))
            count += 1       

        fig, axs = plt.subplots()

        count = 0
        for key in organismList:
            x_vals = []
            y_vals = []
            for i in range(0, len(organismList[key])):
                x_vals.append(float(organismList[key][i][0]))
                y_vals.append(float(organismList[key][i][1]))

            scatter = axs.scatter(x_vals, y_vals, s = 80, label = key, marker=markers[count])

            count += 1

        # produce a legend with the unique colors from the scatter
        legend1 = axs.legend(loc="upper right", title="Off Target Organism")
        legend1 = axs.legend( prop={'size': 6})
        axs.add_artist(legend1)
        axs.set_ylabel('Organism Distance')
        axs.set_xlabel('Off-Target Score')
        axs.set_title('All gRNAs')
        plt.tight_layout()
        plt.show()


if __name__ == "__main__":
   main(sys.argv[1:])