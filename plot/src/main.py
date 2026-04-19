import matplotlib.pyplot as plt

if __name__ == '__main__':
    m_values = []
    s_values = []
    y_values = []
    cos_th_values = []
    weight_values = []

    with open('/home/julian/documents/uni/master/master_thesis/learning/simple_drell_yan/data/events/events.csv') as file:
        for line in file.readlines():
            cols = line.split(", ")
            m_values.append(float(cols[0]))
            s_values.append(float(cols[1]))
            y_values.append(float(cols[2]))
            cos_th_values.append(float(cols[3]))
            weight_values.append(float(cols[4]))

    weight_sum = sum(weight_values)
    weight_values = [w/weight_sum for w in weight_values]
    
    plt.figure()
    plt.hist(m_values, weights=weight_values)
    plt.savefig('/home/julian/documents/uni/master/master_thesis/learning/simple_drell_yan/plot/m_hist.pdf')

    plt.figure()
    plt.hist(s_values, weights=weight_values)
    plt.savefig('/home/julian/documents/uni/master/master_thesis/learning/simple_drell_yan/plot/s_hist.pdf')

    plt.figure()
    plt.hist(y_values, weights=weight_values)
    plt.savefig('/home/julian/documents/uni/master/master_thesis/learning/simple_drell_yan/plot/y_hist.pdf')

    bin_count = 20
    plt.figure()
    plt.hist(cos_th_values, weights=weight_values, bins=[-1 + 2 / (bin_count + 1) * i for i in range(bin_count + 1)])
    plt.savefig('/home/julian/documents/uni/master/master_thesis/learning/simple_drell_yan/plot/cos_th_hist.pdf')
