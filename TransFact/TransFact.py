import numpy as np, click, pandas as pd
from scipy.optimize import minimize
from collections import defaultdict


# Objective function for optimization
def objective_function(params, transmission_matrices, num_targets, num_sources, regularization_strength, valid_mask):
    num_pathogens = len(transmission_matrices)
    n = num_targets * num_sources

    # Extract human vector from params
    human_vector = params[:n].reshape(num_targets, num_sources)

    loss = 0
    prediction = np.zeros((num_pathogens, num_targets, num_sources))
    country_weights = params[n:].reshape([num_pathogens, num_sources])

    # Apply mask: zero out entries for invalid sources
    country_weights *= valid_mask
    country_weights /= (country_weights.sum(1, keepdims=True) + 1e-10)

    for p in range(num_pathogens):
        predicted = human_vector * country_weights[p]
        predicted[transmission_matrices[p] == 0] = 0
        predicted = (predicted.T / (predicted.sum(1) + 1e-10)).T
        prediction[p] = predicted

    loss += np.sum((prediction - transmission_matrices)**2)
    reg_loss = regularization_strength * (np.sum(params[n:] * valid_mask.ravel()) + np.sum(human_vector))

    print("Loss:", loss, "Reg Loss:", reg_loss)
    return loss + reg_loss


def initialize_parameters(transmission_matrices, valid_mask):
    num_pathogens, num_targets, num_sources = transmission_matrices.shape
    n = num_targets * num_sources

    avg_matrix = np.mean(transmission_matrices, axis=0)
    initial_human_vector = avg_matrix.ravel()

    initial_country_weights = []
    for p in range(num_pathogens):
        weights = np.mean(transmission_matrices[p], axis=0) + 1e-3 * np.random.rand(num_sources)
        weights *= valid_mask[p]
        if np.sum(weights) > 0:
            weights = weights / np.max(weights)
        initial_country_weights.append(weights)

    return np.concatenate([initial_human_vector] + initial_country_weights)


# Main function to estimate transmission patterns
def estimate_transmission_patterns(transmission_matrices):
    num_pathogens, num_targets, num_sources = transmission_matrices.shape
    n = num_targets * num_sources

    # Create a mask of valid pathogen-source combinations
    valid_mask = (transmission_matrices.sum(axis=1) > 0).astype(float)

    initial_params = initialize_parameters(transmission_matrices, valid_mask)

    print("Local refinement...")
    result = minimize(
        objective_function,
        initial_params,
        args=(transmission_matrices, num_targets, num_sources, 1e-3, valid_mask),
        method='L-BFGS-B',
        bounds=[(1e-3, 1)] * n + [
            (0 if valid_mask[p, s] == 0 else 1e-3, 1)
            for p in range(num_pathogens)
            for s in range(num_sources)
        ],
        options={
            'maxiter': 500000,
            'disp': True,
            'maxcor': 2000,
            'ftol': 1e-10,
            'gtol': 1e-8
        },
    )

    optimized_params = result.x
    human_vector = optimized_params[:n].reshape(num_targets, num_sources)
    human_vector = (human_vector.T / (human_vector.sum(1) + 1e-10)).T

    country_weights = optimized_params[n:].reshape([num_pathogens, num_sources])
    country_weights *= valid_mask
    country_weights /= (country_weights.sum(0, keepdims=True) + 1e-10)
    return human_vector, country_weights


def parse_transmission(tfile, minimum_involvement) :
    involvements = defaultdict(int)
    transmissions = defaultdict(int)
    with open(tfile, 'rt') as fin :
        for line in fin :
            p = line.strip().split('\t')
            if p[0].startswith('China'): p[0] = 'China'
            if p[1].startswith('China'): p[1] = 'China'
            if p[0] == p[1]: continue
            involvements[p[0]] += 1
            involvements[p[1]] += 1
            transmissions[(p[1], p[0])] += 1

    involvements = {p:v for p, v in involvements.items() if v >= minimum_involvement}
    tmp = defaultdict(int)
    totals = defaultdict(int)
    for (p2, p1), v in transmissions.items() :
        if p1 in involvements and p2 in involvements :
            totals[p2] += v
            tmp[(p2, p1)] += v
    
    transmissions = {}
    involvements = sorted(involvements.keys())
    for i, p1 in enumerate(involvements) :
        for p2 in involvements[:i] :
            if totals[p1] >= minimum_involvement :
                transmissions[(p1, p2)] = tmp.get((p1, p2), 0.1)
            if totals[p2] >= minimum_involvement :
                transmissions[(p2, p1)] = tmp.get((p2, p1), 0.1)
    return transmissions


@click.command()
@click.option('-p', '--prefix', help='prefix for the outputs.', default='transmission', show_default=True)
@click.option('-m', '--minimum_involvement', default=10, show_default=True)
@click.argument('transmission_files', nargs=-1)
def main(prefix, minimum_involvement, transmission_files) :
    transmissions = []
    for tfile in transmission_files :
        transmissions.append([tfile, parse_transmission(tfile, minimum_involvement)])
    
    sources, targets = defaultdict(int), defaultdict(int)
    for tfile, transmission in transmissions :
        s0, t0 = {}, {}
        for (tgt, src), cnt in transmission.items() :
            s0[src], t0[tgt] = 1, 1
        for src in s0: sources[src] += 1
        for tgt in t0: targets[tgt] += 1

    source0 = {}
    for src, cnt in sorted(sources.items()) :
        if cnt >= 2 :
            source0[src] = len(source0)
    target0 = {}
    for tgt, cnt in sorted(targets.items()) :
        if cnt >= 2 :
            target0[tgt] = len(target0)
    sources, targets = source0, target0
    
    transmission_matrices = np.zeros([len(transmissions), len(targets), len(sources)], dtype=float)
    for i, (tfile, transmission) in enumerate(transmissions) :
        for (tgt, src), cnt in transmission.items() :
            if tgt in targets and src in sources :
                transmission_matrices[i, targets[tgt], sources[src]] = cnt
    transmission_matrices /= (np.sum(transmission_matrices, 2).reshape([transmission_matrices.shape[0], transmission_matrices.shape[1], 1])+1e-10)
    
    print("Estimating transmission patterns...")
    V_est, W_list = estimate_transmission_patterns(transmission_matrices)

    vector_migrations = pd.DataFrame(V_est, index=sorted(targets.keys()), columns=sorted(sources.keys()))
    vector_migrations.to_csv(f'{prefix}.vectors.csv')
    
    patho_feats = pd.DataFrame(W_list, index=transmission_files, columns=sorted(sources.keys()))
    patho_feats.to_csv(f'{prefix}.pathogens.csv')
    return



# Main execution
if __name__ == "__main__":
    main()
