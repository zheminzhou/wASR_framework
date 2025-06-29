import click, numpy as np, json, re, collections

treetime = '/home/zhemin/.local/bin/treetime'
wASR = '/titan/softwares/EPHI/wASR/wASR.py'

@click.command()
@click.option('-r', '--result_file')
def main(result_file) :
    results = collections.defaultdict(list)
    with open(result_file, 'rt') as fin :
        for line in fin :
            fname, json_str = line.strip().split(' ', 1)
            res = json.loads(json_str.replace("'", '"'))
            m_rate = 'H' if 'simulation_final0' in fname else 'L'
            a_freq = re.findall(r'x(\d+)_', fname)[0]
            s_freq = re.findall(r'_1_(\d+).ite', fname)[0]
            
            category = f'{m_rate}_{a_freq}_{s_freq}'
            for tool, r in res.items() :
                results[(category, tool)].append((r[0]+r[2])/(r[0]+r[1]+r[2]+r[3]))
                # results[(category, tool)].append(r[0]/(r[0]+r[1]+0.5))
    for (category, tool), dat in results.items() :
        print(f'{category}\t{tool}\t{np.mean(dat)}\t{np.std(dat)}')


if __name__ == '__main__' :
    main()
