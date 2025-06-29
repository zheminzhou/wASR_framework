import click, numpy as np, ete3_extensions, os, subprocess, re, shutil, pandas as pd

treetime = '/home/zhemin/miniconda3/envs/p312/bin/treetime'
wASR = '/titan/softwares/EPHI/wASR/wASR.py'

ape_weighted_script = open("wACE.R").read() + '''
tree <- read.tree("tree.nwk")
tree$node.label[1] <- "NODE_0000000"
states_data <- read.csv("states.csv")

states <- states_data$state
names(states) <- states_data$ID
states <- as.factor(states[tree$tip.label])

states_weight <- read.csv("states.weight")
weights <- states_weight$weight
names(weights) <- states_weight$state
weights <- weights[states]
weights <- weights/mean(weights)
# weights <- as.factor(weights)

ace_result <- ace_weighted(states, tree, type = "discrete", model = "SYM", weights=weights)
node_states <- apply(ace_result$lik.anc, 1, function(x) colnames(ace_result$lik.anc)[which.max(x)])
write.csv(data.frame(
  node = (length(tree$tip.label) + 1):(length(tree$tip.label) + tree$Nnode),
  label = tree$node.label,
  most_likely_state = node_states,
  probability = apply(ace_result$lik.anc, 1, max)
), "SYM_result.csv", row.names = FALSE)
'''

ape_script = '''library(ape)
tree <- read.tree("tree.nwk")
tree$node.label[1] <- "NODE_0000000"
states_data <- read.csv("states.csv")

states <- states_data$state
names(states) <- states_data$ID
states <- as.factor(states[tree$tip.label])

ace_result <- ace(states, tree, type = "discrete", model = "SYM")
node_states <- apply(ace_result$lik.anc, 1, function(x) colnames(ace_result$lik.anc)[which.max(x)])
write.csv(data.frame(
  node = (length(tree$tip.label) + 1):(length(tree$tip.label) + tree$Nnode),
  label = tree$node.label,
  most_likely_state = node_states,
  probability = apply(ace_result$lik.anc, 1, max)
), "SYM_result.csv", row.names = FALSE)

# ace_result <- ace(states, tree, type = "discrete", model = "ARD")
# node_states <- apply(ace_result$lik.anc, 1, function(x) colnames(ace_result$lik.anc)[which.max(x)])
# write.csv(data.frame(
#   node = (length(tree$tip.label) + 1):(length(tree$tip.label) + tree$Nnode),
#   label = tree$node.label,
#   most_likely_state = node_states,
#   probability = apply(ace_result$lik.anc, 1, max)
# ), "ARD_result.csv", row.names = FALSE)
'''


def eval_treetime(confidence, ground_truth) :
    dat = pd.read_csv(confidence)
    res = [0, 0, 0, 0]
    for d in dat.values :
        if d[0] in ground_truth :
            trait, prob = ('x0', d[1]) if d[1] >= 0.5 else ('x1', d[2])
            i = 0 if prob >= 0.7 else 2
            i += 0 if trait == ground_truth[d[0]] else 1
            res[i] += 1
    return res


@click.command()
@click.option('-n', '--nexus')
@click.option('-t', '--task', default='all')
def main(nexus, task) :
    tasks = task.split(',')
    results = {}
    # prepare folder
    odir = nexus.rsplit('.', 1)[0]
    if not os.path.isdir(odir) :
        os.makedirs(odir)
    
    # set up required files
    simtree = ete3_extensions.read_nexus(nexus)[0]
    simtree.name = 'NODE_0000000'
    ground_truth = {node.name:node.annotations['states'] for node in simtree.traverse() if not node.is_leaf()}
    traits = []
    for node in simtree.traverse() :
        node.dist *= 0.1
    
    for tip in simtree.get_leaves() :
        traits.append([tip.name, tip.annotations['states']])
    
    with open(os.path.join(odir, 'tree.nwk'), 'wt') as fout :
        fout.write(simtree.write(format=1))
    
    with open(os.path.join(odir, 'states.csv'), 'wt') as fout :
        fout.write('ID,state\n')
        for trait in traits :
            fout.write(f'{trait[0]},{trait[1]}\n')
    
    freq = [float(f) for f in re.findall(r'_x(\d+)_(\d+)_sim', nexus)[0]]

    # run treetime
    if 'all' in tasks or 'treetime' in tasks :
        subprocess.run(f'{treetime} mugration --states states.csv --tree tree.nwk --confidence --out treetime.out'.split(), cwd=odir, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        res = eval_treetime(f'{odir}/treetime.out/confidence.csv', ground_truth)
        results['treetime'] = res

    # run wASR
    if 'all' in tasks or 'wASR' in tasks :
        state_count = np.unique([t[1] for t in traits], return_counts=True)[1]
        base_num = 0 # np.log(len(traits))
        weights = [freq[0]/(state_count[0]+base_num), freq[1]/(state_count[1]+base_num)]
        #print(weights, weights[0]/weights[1])
        with open(os.path.join(odir, 'states.weight'), 'wt') as fout :
            fout.write('state,weight\n')
            for i, w in enumerate(weights) :
                fout.write(f'x{i},{w}\n')

        subprocess.run(f'python {wASR} -w states.weight --states states.csv --tree tree.nwk --confidence --out wASR.out'.split(), cwd=odir, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        res = eval_treetime(f'{odir}/wASR.out/confidence.csv', ground_truth)
        results['wASR'] = res

    # run ape
    if 'all' in tasks or 'ape' in tasks :
        with open(os.path.join(odir, 'ape.script'), 'wt') as fout :
            fout.write(ape_script)
        subprocess.run('/home/zhemin/miniconda3/envs/R/bin/R -f ape.script'.split(), cwd=odir, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        
        res = eval_ape(os.path.join(odir, 'SYM_result.csv'), ground_truth)
        results['APE_SYM'] = res
        # res = eval_ape(os.path.join(odir, 'ARD_result.csv'), ground_truth)
        # results['APE_ARD'] = res
    
    # run weighted_ape
    if 'all' in tasks or 'wAPE' in tasks :
        with open(os.path.join(odir, 'ape_weighted.script'), 'wt') as fout :
            fout.write(ape_weighted_script)
        subprocess.run('/home/zhemin/miniconda3/envs/R/bin/R -f ape_weighted.script'.split(), cwd=odir, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        
        res = eval_ape(os.path.join(odir, 'SYM_result.csv'), ground_truth)
        results['wAPE_SYM'] = res
        # res = eval_ape(os.path.join(odir, 'ARD_result.csv'), ground_truth)
        # results['wAPE_ARD'] = res
    
    # outputs
    # shutil.rmtree(odir)
    print(nexus, results)


def eval_ape(res_file, ground_truth) :
    res = pd.read_csv(res_file).values
    d = [0, 0, 0, 0]
    for r in res :
        i = 0 if r[3] >= 0.7 else 2
        i += 1 if r[2] != ground_truth[r[1]] else 0
        d[i] += 1
    return d



if __name__ == '__main__' :
    main()
