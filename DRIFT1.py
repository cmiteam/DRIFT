import os
import csv
import math
import random
import numpy as np
from PIL import Image
import scipy.stats as stats
from bitarray import bitarray
import matplotlib.pyplot as plt
from collections import defaultdict

random.seed()
data_directory = "data"
results_directory = "results"

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

def run_population_model(model):

    free_params = setup_free_params(model)
    chromosome_arm_data, free_params['numbits'] = load_chromosome_data(model["multiplier"])
    death_risk = load_actuarial_table()
    model_name = model["model_id"]

    print("Model parameters: ", model)
    print()
    for run in range(1, model["num_runs"] + 1):
        print ("Run ", run)

        IndData = {}
        chromosomes = {}
        mutations = {}
        mutation_hist = {i: 0 for i in range(-1000, 1001)}
        setup_output_files(model, run)
        initialize_population(IndData, model, free_params, chromosomes, mutations)
        things_to_plot_during_run = setup_plot(model)
        fig, ax1 = plt.subplots()
        ax2 = ax1.twinx()
        tracking = {'marriages': 0, 'births': 0, 'random_deaths': 0, 'cull_deaths': 0}

        for year in range(0, model["end_year"] + 1):

            if year == model["seed_year"]:
                Innoculate_random_person(IndData, chromosomes, free_params)
            babymommas, births = putemintheoven(IndData, year, model)
            tracking['births'] += births
            mutation_hist = birth(IndData, model, free_params, babymommas, year, chromosomes, chromosome_arm_data, mutations, mutation_hist, free_params['numbits'])
            manlist, womanlist = list_availables(model, IndData, year)
            tracking['marriages'] += setup_marriages(manlist, womanlist, IndData, model)
            lifesucks, culled = BumpPeopleOff(IndData, chromosomes, mutations, model, free_params, death_risk, year, run)
            tracking['random_deaths'] += lifesucks
            tracking['cull_deaths'] += culled

            if len(IndData) == 0:
                print ("Population Extinct!")
                break
                # TO DO: save population data midrun if population goes extinct

            if year % model['save_interval'] == 0:
                Save(run, year, IndData, chromosomes, model, free_params, tracking, things_to_plot_during_run, fig, ax1, ax2, mutations)
                tracking = {'marriages': 0, 'births': 0, 'random_deaths': 0, 'cull_deaths': 0}
                if model['track_DNA'] and model['every_genome_map']:
                    save_population_genome_map(f"{model_name}-{year}", free_params['numbits'], IndData, chromosomes, chromosome_arm_data)
            free_params['lastpopsize'] = len(IndData)

        if model['track_mutations']:
            if model['mutation_hist']:
                save_mutation_histogram(mutations, mutation_hist, model_name, run)
            if model['mutation_map']:
                save_population_mutation_map(model_name, free_params['numbits'], IndData, mutations, chromosome_arm_data)

        if model['track_DNA'] and model['genome_map']:
            save_population_genome_map(model_name, free_params['numbits'], IndData, chromosomes, chromosome_arm_data)

        if model['track_dead']:
            save_still_living_people(IndData, model['model_id'], model['end_year'], model['num_runs'])

    print ('Done')

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

def Save(run, year, IndData, chromosomes, model, free_params, tracking, things_to_plot_during_run, fig, ax1, ax2, mutations):

    model_name = model["model_id"]
    filename = os.path.join(results_directory, f"{model_name}-{run} results.csv")

    numinds = len(IndData)

    Y_descends, mt_descends, genealo_descends, genetic_descends, num_blocks, num_centromeres = 0, 0, 0, 0, 0, 0
    tot_blocks, av_block_size, sd_block_size = 0, 0, 0
    perc_seed_genome_retained, av_seed_genome_coverage, av_heterozygosity = 0, 0, 0
    av_ind_fitness, av_bin_fitness, num_mutations, av_mutations_per_ind, av_mutations_per_bin = 0, 0, 0, 0, 0

    if model["track_DNA"] == 1:
        Y_descends, mt_descends, genealo_descends, genetic_descends, num_blocks, num_centromeres = calculate_misc_stats(IndData)
        tot_blocks, av_block_size, sd_block_size = calculate_block_stats(chromosomes)
        perc_seed_genome_retained, av_seed_genome_coverage, av_heterozygosity = calculate_genetic_stats(chromosomes, free_params["numbits"], numinds)

    if model["track_mutations"] == 1:
        av_ind_fitness, av_bin_fitness, num_mutations, av_mutations_per_ind, av_mutations_per_bin = calculate_fitness_stats(IndData, free_params["numbits"])
        av_ind_fitness = format(av_ind_fitness, '.01f')
        av_bin_fitness = format(av_bin_fitness, '.01f')
        av_mutations_per_ind = format(av_mutations_per_ind, '.1f')
        av_mutations_per_bin = format(av_mutations_per_bin, '.1f')

    # Save current model status
    with open(filename, mode='a', newline='') as results_file:
        results_writer = csv.writer(results_file)
        results_writer.writerow([run, year, numinds, tracking['marriages'], tracking['births'], tracking['random_deaths'], tracking['cull_deaths'], genetic_descends, genealo_descends, Y_descends, mt_descends, num_centromeres, tot_blocks, av_block_size, sd_block_size, av_ind_fitness, av_bin_fitness, num_mutations, av_mutations_per_ind, av_mutations_per_bin, perc_seed_genome_retained, av_seed_genome_coverage, av_heterozygosity])

    # Track progress on screen
    numinds = len(IndData)
    print(f"Y:{year} N:{numinds} m:{tracking['marriages']} b:{tracking['births']} r:{tracking['random_deaths']} c:{tracking['cull_deaths']} ge:{genetic_descends} go:{genealo_descends} Y:{Y_descends} mt:{mt_descends} cs:{num_centromeres} mI:{free_params['indID']} bl:{num_blocks} abl:{av_block_size} sbl:{sd_block_size} FI:{av_ind_fitness} FB:{av_bin_fitness} NM:{num_mutations} MI:{av_mutations_per_ind} MB:{av_mutations_per_bin}")

    # Update onscreen graph
#    marriages = tracking['marriages']
#    births = tracking['births']
#    random_deaths = tracking['random_deaths']
#    cull_deaths = tracking['cull_deaths']
#    max_ID = free_params['indID']
    for variable in model["plot"]:
        value = locals()[variable]
        things_to_plot_during_run[variable].append(value)
    things_to_plot_during_run['year'].append(year)

    variables_to_skip = ['perc_seed_genome_retained', 'av_seed_genome_coverage', 'av_heterozygosity', 'av_fitness_per_ind', 'av_fitness_per_bin']     # these are generally not larger than 1.0, so they go on the second y axis
    ax1.clear()
    ax2.clear()

    for variable in model["plot"]:
        if variable != 'year' and variable not in variables_to_skip:
            ax1.plot(things_to_plot_during_run["year"], things_to_plot_during_run[variable], label=variable)
    for variable in variables_to_skip:
        if variable in things_to_plot_during_run:
            ax2.plot(things_to_plot_during_run["year"], things_to_plot_during_run[variable], label=variable, linestyle='dashed')
            ax2.legend(loc='upper right')
    ax1.set_xlabel('Year')
    ax1.legend(loc='upper left')
    plt.draw()
    plt.pause(0.1)

    plot_geography = 0
    if plot_geography == 1:
        coordinates = [(ind['lat'], ind['lon']) for ind in IndData.values()]
        latitudes, longitudes = zip(*coordinates)
        plt.scatter(latitudes, longitudes, label='Individuals')
        circle = plt.Circle((0, 0), 0.5, edgecolor='r', facecolor='none', linestyle='dashed')
        plt.gca().add_patch(circle)
        plt.xlim(-0.6, 0.6)
        plt.ylim(-0.6, 0.6)
        plt.xlabel('Latitude')
        plt.ylabel('Longitude')
        plt.title('Geographical Locations of Individuals')
        plt.legend()
        plt.draw()
        plt.pause(0.1)
    
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

def save_still_living_people (IndData, model_name, year, run):

    dead_people_data = ''
    for ind in sorted(IndData.keys()):
        dead_people_data += dead_string(ind, IndData, year)
    filename = os.path.join(results_directory, f"{model_name}-{run} deaths.csv")
    with open(filename, mode='a', newline='') as tracked_dead_file:
        tracked_dead_file.write(dead_people_data)

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

def save_mutation_histogram(mutations, mutation_hist, model_name, run):

    circulating_mutations = {i: 0 for i in range(-1000, 1001)}
    for ind in mutations:
        for copy in range (2):
            for position in range(len(mutations[ind][copy])):
                fitness_list = mutations[ind][copy][position]
                for fitness_effect in fitness_list:
                    circulating_mutations[fitness_effect] += 1

    bins = np.arange(-1000, 1000 + 1, 1)
    filename = os.path.join(results_directory, f"{model_name}-{run} mutation_histogram.csv")
    with open(filename, mode='w', newline='') as results_file:
        writer = csv.writer(results_file)
        writer.writerow(['Bin', 'All_mutations', 'Circulating_mutations'])
        for bin in bins:
            adjusted_bin = bin / 1000
            writer.writerow([adjusted_bin, mutation_hist[bin], circulating_mutations[bin]])

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

def save_population_genome_map(model_name, numbits, IndData, chromosomes, chromosome_arm_data):

    image_width = numbits
    image_height = len(IndData) * 2 + 10
    image = Image.new('RGB', (image_width, image_height), 'white')
    pixels = image.load()

    for x in range(numbits):
        for j in range(10):
            pixels[x, j] = (255, 255, 255)

    for chrom in range(1, len(chromosome_arm_data) + 1):
#        pstart, plen = chromosome_arm_data[chrom][0]
#        qstart, qlen = chromosome_arm_data[chrom][1]
        pstart = chromosome_arm_data[chrom]['p']['start']
        qstart = chromosome_arm_data[chrom]['q']['start']
        qlen = chromosome_arm_data[chrom]['q']['length']
        # draw telomere
        for x in range(pstart, pstart + 3):
            for j in range(10):
               pixels[x, j] = (0, 255, 0)
        # draw centromere
        for x in range(qstart - 3, qstart + 3):
            for j in range(10):
                pixels[x, j] = (0, 0, 255)
        # draw telomere
        for x in range(qstart + qlen - 4, qstart + qlen - 1):
            for j in range(10):
                pixels[x, j] = (0, 255, 0)
        # draw separator
        for j in range(10):
            pixels[qstart + qlen - 1, j] = (0, 0, 0)

    row = 10
    for ind in IndData:
        if ind in chromosomes:
            for copy in range (2):
                chromosome = chromosomes[ind][copy]
                for j in range(numbits):
                    if chromosome[j]:
                        pixels[j, row] = (255,0,0)
                    else:
                        pixels[j, row] = (0,0,0)
                row += 1
        else:
            for copy in range (2):
                for j in range(numbits):
                   pixels[j, row] = (0,0,0)
                row += 1

    new_image_height = row
    resized_image = image.resize((image_width, new_image_height))
    filename = os.path.join(results_directory, f"{model_name} genome_map.png")
    resized_image.save(filename, format='PNG')

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

def save_population_mutation_map(model_name, image_width, IndData, mutations, chromosome_arm_data):

    image_height = len(IndData) * 2 + 10
    image = Image.new('RGB', (image_width, image_height), 'white')
    pixels = image.load()
    for chrom in range(1, len(chromosome_arm_data) + 1):
#        pstart, plen = chromosome_arm_data[chrom][0]
#        qstart, qlen = chromosome_arm_data[chrom][1]
        pstart = chromosome_arm_data[chrom]['p']['start']
        qstart = chromosome_arm_data[chrom]['q']['start']
        qlen = chromosome_arm_data[chrom]['q']['length']
        # draw telomere
        for x in range(pstart, pstart + 3):
            for j in range(10):
               pixels[x, j] = (0, 255, 0)
        # draw centromere
        for x in range(qstart - 3, qstart + 3):
            for j in range(10):
                pixels[x, j] = (0, 0, 255)
        # draw telomere
        for x in range(qstart + qlen - 4, qstart + qlen - 1):
            for j in range(10):
                pixels[x, j] = (0, 255, 0)
        # draw separator
        for j in range(10):
            pixels[qstart + qlen - 1, j] = (0, 0, 0)

    row = 10
    for ind in IndData:
        if ind in mutations:
            for copy in range (2):
                for j in range(image_width):
                    fitness, mutation_count = get_bin_fitness(ind, copy, j, mutations)
                    if mutation_count == 0:
                        pixels[j, row] = (0,0,0)
                    else:
                        if fitness == 1:
                            pixels[j, row] = (255,255,255)
                        if fitness > 1:
                            pixels[j, row] = (0,255,0)
                        if fitness < 1:
                            pixels[j, row] = (255,0,0)
                row += 1
        else:
            for copy in range (2):
                for j in range(image_width):
                   pixels[j, row] = (0,0,0)
                row += 1

    filename = os.path.join(results_directory, f"{model_name} mutation_map.png")
    image.save(filename, format='PNG')

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

def get_bin_fitness(ind, copy, position, mutations):
    fitness_list = mutations[ind][copy][position]
    fitness = 1
    mutation_count = len(fitness_list)
    for j in range(mutation_count):
        fitness += fitness_list[j] / 1000
    return fitness, mutation_count

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

def list_availables(model, IndData, year):

    grooms = []
    brides = []
    for ind in IndData:
        married = IndData[ind]["marriage_state"]
        if married == -1:
            age = year - IndData[ind]["birth_year"]
            if age > model["maturity"]:
                sex = IndData[ind]["sex"]
                lifespan  = IndData[ind]["lifespan"]
                if sex == 0:
                    grooms.append(ind)
                elif age < lifespan * model["menopause"]:    # Old ladies don't remarry
                    brides.append(ind)

    return grooms, brides

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

def putemintheoven(IndData, year, model):

    pregnant_couples = []
    for ind in IndData:
        sex      = IndData[ind]["sex"]
        married  = IndData[ind]["marriage_state"]
        age      = year - IndData[ind]["birth_year"]
        lifespan = IndData[ind]["lifespan"]
        if sex == 1 and married > -1 and age < lifespan * model["menopause"]:
            if IndData[ind]["year_of_last_birth"] + model["spacing"] <= year:
                if random.randint(0, model["birth_prob"] - 1) == 0:
                    mom = ind
                    dad = married
                    fitness = 1
                    if model["track_mutations"] == 1 and model["selection"] == "birth":
                        fitness = (IndData[dad]["fitness"] + IndData[mom]["fitness"]) / 2
                    chance = random.random()
                    if chance < fitness:
                        b = {'dad': dad, 'mom': mom}
                        pregnant_couples.append(b)

    return pregnant_couples, len(pregnant_couples)

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

def create_child(year, IndData, dad, mom, child, model):

    life = int((IndData[dad]['lifespan'] + IndData[mom]['lifespan']) / 2 * model['lifespan_drop'])
    if life < model['min_lifespan']:
        life = model['min_lifespan']

    IndData[child] = {
        'dad': dad,
        'mom': mom,
        'sex': random.randint(0, 1),
        'birth_year': year,
        'lifespan': life,
        'fitness': 1,
        'marriage_state': -1,
        'lat': round((IndData[dad]['lat'] + IndData[mom]['lat']) / 2, 2),
        'lon': round((IndData[dad]['lon'] + IndData[mom]['lon']) / 2, 2)
    }

    IndData[mom]['year_of_last_birth'] = year
    if 'numbirths' not in IndData[mom]:
        IndData[mom]['numbirths'] = 0
    IndData[mom]['numbirths'] += 1

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

def createmask (sex, free_params, chromosome_arm_data):

    mask  = bitarray(free_params["numbits"])
    cents = bitarray(len(chromosome_arm_data) + 1)
    mask.setall(0)
    cents.setall(0)
    for chrom in range(1, len(chromosome_arm_data) + 1):
        pstart = chromosome_arm_data[chrom]['p']['start']
        qstart = chromosome_arm_data[chrom]['q']['start']
        plen = chromosome_arm_data[chrom]['p']['length']
        qlen = chromosome_arm_data[chrom]['q']['length']
        ploc = random.randint(0, plen - 1)
        qloc = random.randint(0, qlen - 1)
        which_copy = random.randint(0, 1)
        if which_copy:
            mask[pstart + ploc:qstart + qloc] = True
        else:
            mask[pstart:pstart + ploc] = True
            mask[qstart + qloc:qstart + qlen] = True
            cents[chrom] = 1

    if sex == 0:
        xstart = chromosome_arm_data[23]['p']['start']
        xend =  chromosome_arm_data[23]['q']['start'] +  chromosome_arm_data[23]['q']['length']
        mask[xstart:xend] = False    

    return mask, cents

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

def meiosis(parent, child, chromosomes, mask, copy):
    parent_copy_0 = bitarray(len(mask))
    parent_copy_1 = bitarray(len(mask))
    child_copy = bitarray(len(mask))
    parent_copy_0.setall(0)
    parent_copy_1.setall(0)
    child_copy.setall(0)
    if len(chromosomes[parent][0]) == len(mask):
        parent_copy_0 = chromosomes[parent][0]
    if len(chromosomes[parent][1]) == len(mask):
        parent_copy_1 = chromosomes[parent][1]
    child_copy = (mask & parent_copy_0) | (~mask & parent_copy_1)
    chromosomes[child][copy] = child_copy

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

def count_alleles (ind, chromosomes):

    allele_count = 0
    allele_count = chromosomes[ind][0].count(True)
    allele_count += chromosomes[ind][1].count(True)
    return allele_count

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

def inherit_centromeres (cents1, cents2, IndData, dad, mom, child, chromosome_arm_data):

    def update_child_centromeres(IndData, child, parent, chromosome_arm_data, offset, cents):
        for chrom in range(1, len(chromosome_arm_data) + 1):
            if cents[chrom] == 0:
                IndData[child]['centromeres'][chrom * 2 + offset] = IndData[parent]['centromeres'][chrom * 2]
            else:
                IndData[child]['centromeres'][chrom * 2 + offset] = IndData[parent]['centromeres'][chrom * 2 + 1]

    if 'centromeres' in IndData[dad] or 'centromeres' in IndData[mom]:
        IndData[child]['centromeres'] = bitarray(48)
        IndData[child]['centromeres'].setall(0)
        if 'centromeres' in IndData[dad]:
            update_child_centromeres(IndData, child, dad, chromosome_arm_data, 0, cents1)
        if 'centromeres' in IndData[mom]:
            update_child_centromeres(IndData, child, mom, chromosome_arm_data, 1, cents2)

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

def inherit_mutations(mask1, mask2, dad, mom, child, mutations, numbits):

    mutations[child][0] = [mutations[dad][int(mask1[i])][i] for i in range(numbits)]
    mutations[child][1] = [mutations[mom][int(mask2[i])][i] for i in range(numbits)]
    return mutations

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

def generate_new_mutations(ind, model, free_params, mutations, mutation_hist):

    def get_mutation_effect(model):
        mutation_effect = 0
        is_mutation_non_neutral = random.random()
        if is_mutation_non_neutral >= model['f_neutral']:
            mutation_effect = stats.weibull_min.rvs(model['shape'], scale=model['scale'])
            mutation_effect = mutation_effect / model['Weibull_adj']
            is_mutation_deleterious = random.random()
            if is_mutation_deleterious > model['f_beneficial']:
                mutation_effect = -mutation_effect
        return mutation_effect

    numbits = free_params['numbits']
    num_new_mutations = np.random.poisson(model['mu'])
    for _ in range(num_new_mutations):
        position = random.randint(0, numbits - 1)
        mutation_effect = get_mutation_effect(model)
        free_params['mutID'] += 1
        integer_value = int(mutation_effect * 1000)
        mutation_hist[integer_value] += 1
        which_copy = random.randint(0, 1)
        mutations[ind][which_copy][position].append(integer_value)
    return mutations, mutation_hist
    
# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

def count_fitness_and_mutations(mutations, ind):

    mutation_count = 0
    fitness = 1
    for position in range(len(mutations[ind][0])):
        mutation_count += len(mutations[ind][0][position])
        if mutation_count > 0:
            fitness += sum(mutations[ind][0][position]) / (1000 * mutation_count)
    for position in range(len(mutations[ind][1])):
        mutation_count += len(mutations[ind][1][position])
        if mutation_count > 0:
            fitness += sum(mutations[ind][0][position]) / (1000 * mutation_count)

    return mutation_count, fitness

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

def birth(IndData, model, free_params, birthlist, year, chromosomes, chromosome_arm_data, mutations, mutation_hist, numbits):

    for birth_info in birthlist:
        dad = birth_info['dad']
        mom = birth_info['mom']
        free_params['indID'] += 1
        child = free_params['indID']
        create_child(year, IndData, dad, mom, child, model)

        if model['track_DNA'] or model['track_mutations']:
            mask1, cents1 = createmask(0, free_params, chromosome_arm_data)
            mask2, cents2 = createmask(1, free_params, chromosome_arm_data)

        if model['track_DNA']:
            chromosomes[child] = {}
            chromosomes[child][0] = bitarray(numbits)
            chromosomes[child][1] = bitarray(numbits)
            chromosomes[child][0].setall(0)
            chromosomes[child][1].setall(0)
            if IndData[dad]['allele_count'] > 0:
                meiosis(dad, child, chromosomes, mask1, 0)
            if IndData[mom]['allele_count'] > 0:
                meiosis(mom, child, chromosomes, mask2, 1)
            IndData[child]['allele_count'] = count_alleles(child, chromosomes)
            if IndData[child]['allele_count'] < 1 or IndData[child]['sex'] == 1:
                del chromosomes[child]
            else:
                IndData[child]['num_blocks'] = count_blocks(chromosomes, child)
            
            inherit_centromeres (cents1, cents2, IndData, dad, mom, child, chromosome_arm_data)

            IndData[child]['Y_gens'] = -1
            if IndData[dad]['Y_gens'] > -1 and IndData[child]['sex'] == 0:
                IndData[child]['Y_gens'] = IndData[dad]['Y_gens'] + 1
            IndData[child]['mt_gens'] = -1
            if IndData[mom]['mt_gens'] > -1:
                IndData[child]['mt_gens'] = IndData[mom]['mt_gens'] + 1
            
            IndData[child]['min_genealo_gens'] = -1
            min_genealo = IndData[dad]['min_genealo_gens']
            if IndData[mom]['min_genealo_gens'] > min_genealo:
                min_genealo = IndData[mom]['min_genealo_gens']
            if min_genealo > -1:
                IndData[child]['min_genealo_gens'] = min_genealo + 1

            IndData[child]['max_genealo_gens'] = -1
            max_genealo = IndData[dad]['max_genealo_gens']
            if IndData[mom]['max_genealo_gens'] > max_genealo:
                max_genealo = IndData[mom]['max_genealo_gens']
            if max_genealo > -1:
                IndData[child]['max_genealo_gens'] = max_genealo + 1

        if model['track_mutations']:
            mutations[child] = [[[[], []] for _ in range(numbits)] for _ in range(2)]
            mutations = inherit_mutations(mask1, mask2, dad, mom, child, mutations, numbits)
            mutations, mutation_hist = generate_new_mutations(child, model, free_params, mutations, mutation_hist)
            IndData[child]['mutations'], IndData[child]['fitness'] = count_fitness_and_mutations(mutations, child)
    return mutation_hist

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

def count_blocks(chromosomes, ind):

    num_blocks = 0
    for copy in range(2):
        bit_array = chromosomes[ind][copy]
        bit_string = bit_array.to01()
        if bit_string == '0' * len(bit_string):  # if the bit array contains only unset bits
            continue
        substrings = bit_string.split("0")
        substrings = [substring for substring in substrings if substring]
        num_blocks += len(substrings)

    return num_blocks

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

def setup_marriages(manlist, womanlist, IndData, model):

    newmarriages = 0
    if model['random_mating'] == 1:
        newmarriages = random_mariages(manlist, womanlist, IndData, model)
    else:
        newmarriages = non_random_mariages(manlist, womanlist, IndData, model)
    return newmarriages

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

def random_mariages(manlist, womanlist, IndData, model):

    random.shuffle(manlist)
    random.shuffle(womanlist)
    men = len(manlist)
    women = len(womanlist)
    max_count = min(men, women)

    for i in range(max_count):
        man = manlist[i]
        woman = womanlist[i]
        IndData[man]['marriage_state'] = woman
        IndData[woman]['marriage_state'] = man
        IndData[woman]['year_of_last_birth'] = -model['spacing']

    return max_count

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

def non_random_mariages(manlist, womanlist, IndData, model):

    # This will pick the man-woman pair closest in age (preferring woman younger than men unless none are avaialable),
    # but marriages will not happen if they are > randommating apart in x,y space

    random.shuffle(manlist)
    random.shuffle(womanlist)
    men = len(manlist)
    women = len(womanlist)
    max_count = min(men, women)

    for i in range(max_count):
        man = manlist[i]
        x1, y1 = IndData[man]['lat'], IndData[man]['lon']

        eligible_females = {}
        for woman in womanlist:
            x2, y2 = IndData[woman]['lat'], IndData[woman]['lon']
            distance = math.sqrt((x2 - x1)**2 + (y2 - y1)**2)
            if distance <= model['random_mating']:
               eligible_females[woman] = IndData[woman]['birth_year'] - IndData[man]['birth_year']

        if eligible_females:
            positive_age_diff = {k: v for k, v in eligible_females.items() if v > 0}
            negative_age_diff = {k: v for k, v in eligible_females.items() if v <= 0}
            if positive_age_diff:
                woman = min(positive_age_diff.items(), key=lambda x: x[1])
                IndData[man]['marriage_state'] = woman
                IndData[woman]['marriage_state'] = man
                IndData[woman]['year_of_last_birth'] = 0
                womanlist.remove(woman)
            else:
                woman  = min(negative_age_diff.items(), key=lambda x: x[1])
                IndData[man]['marriage_state'] = woman
                IndData[woman]['marriage_state'] = man
                IndData[woman]['year_of_last_birth'] = 0
                womanlist.remove(woman)

    return max_count

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

def BumpPeopleOff(IndData, chromosomes, mutations, model, free_params, death_risk, year, run):

    random_deaths, culled_deaths, deaths = 0, 0, 0
    dead_people = []
    dead_people_data = ""

    # Random actuarial deaths
    for ind in IndData:
        age = year - IndData[ind]["birth_year"]
        if ind == free_params["seed"]:
            if age >= IndData[ind]["lifespan"]:
                dead_person_string = dead_string(ind, IndData, year)
                dead_people_data += dead_person_string
                dead_people.append(ind)
            else:
                continue
        fitness = 1
        if model["track_mutations"] == 1 and model["selection"] == "annual":
            fitness = 1 + IndData[ind]["fitness"] / 1000
        die = random.random()
        risk_modification = model["min_lifespan"] / IndData[ind]["lifespan"]
        if 0 < age < 5 and ind != free_params["seed"]:
            if die < 0.0203 * risk_modification + 1 - fitness:
                IndData[ind]['age_at_death'] = age
                IndData[ind]['cause_of_death'] = 'r'
                deaths += 1
                random_deaths += 1
                if model["track_dead"]:
                    dead_person_string = dead_string(ind, IndData, year)
                    dead_people_data += dead_person_string
                dead_people.append(ind)
        else:
            # The actuarial table is in increments of 5, so:
            age_group = int( (age / IndData[ind]["lifespan"]) * model["min_lifespan"] / 5) * 5
            if age_group > 85: # The actuarial table stops at 85, realy old people all have the same probability of dying
                age_group = 85
            if die < death_risk[age_group] * risk_modification + 1 - fitness and ind != free_params["seed"]:
                IndData[ind]['age_at_death'] = age
                IndData[ind]['cause_of_death'] = 'r'
                deaths += 1
                random_deaths += 1
                if model["track_dead"] == 1:
                    dead_person_string = dead_string(ind, IndData, year)
                    dead_people_data += dead_person_string
                dead_people.append(ind)

    RIP(dead_people, IndData, model, chromosomes, mutations)

    max_pop_size = model["max_pop_size"]
    if model['bottleneck_start'] <= year and model['bottleneck_end'] >= year:
        max_pop_size = model['bottleneck_size']

    # trim excess population
    if model['max_breeding_inds'] == 0:
        while len(IndData) > max_pop_size:
            ind = random.choice(list(IndData.keys()))
            if ind == free_params["seed"]:
                continue
            else:
                age = year - IndData[ind]["birth_year"]
                IndData[ind]['age_at_death'] = age
                IndData[ind]['cause_of_death'] = 'c'
                deaths += 1
                culled_deaths += 1
                if model["track_dead"] == 1:
                    dead_people_data += dead_string(ind, IndData, year)
                dead_people = []
                dead_people.append(ind)
                RIP(dead_people, IndData, model, chromosomes, mutations)

    # Tamp down population growth rate
    expected_numinds = free_params["lastpopsize"] * model["max_growth_rate"]
    while len(IndData) > expected_numinds:
        ind = random.choice(list(IndData.keys()))
        if ind == free_params["seed"]:
            continue
        else:
            age = year - IndData[ind]["birth_year"]
            IndData[ind]['age_at_death'] = age
            IndData[ind]['cause_of_death'] = 'c'
            deaths += 1
            culled_deaths += 1
            if model["track_dead"] == 1:
                dead_people_data += dead_string(ind, IndData, year)
            dead_people = []
            dead_people.append(ind)
            RIP(dead_people, IndData, model, chromosomes, mutations)

    # Reduce to specified number of breeding individuals
    if model['max_breeding_inds']:
        breeders = count_breeding_individuals(IndData, year, model)
        while breeders > max_pop_size:
            ind = random.choice(list(IndData.keys()))
            if ind in free_params["seed"]:
                continue
            else:
                age = year - IndData[ind]["birth_year"]
                IndData[ind]['age_at_death'] = age
                IndData[ind]['cause_of_death'] = 'c'
                deaths += 1
                culled_deaths += 1
                if model["track_dead"] == 1:
                    dead_people_data += dead_string(ind, IndData, year)
                dead_people = []
                dead_people.append(ind)
                RIP(dead_people, IndData, model, chromosomes, mutations)
                breeders = count_breeding_individuals(IndData, year, model)

    if model["track_dead"] == 1:
        model = model["model_id"]
        filename = os.path.join(results_directory, f"{model}-{run} deaths.csv")
        with open(filename, mode='a', newline='') as tracked_dead_file:
            tracked_dead_file.write(dead_people_data)
    return random_deaths, culled_deaths

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

def dead_string(ind, IndData, year):

    if ind in IndData:
        ind_info = IndData[ind]
        ind_info['dad'] = ind_info.get('dad', -1)
        ind_info['mom'] = ind_info.get('mom', -1)
        ind_info['numbirths'] = ind_info.get('numbirths', -1)
        ind_info['Ygens'] = ind_info.get('Y_gens', -1)
        ind_info['MTGens'] = ind_info.get('mt_gens', -1)
        ind_info['MaxGenealGens'] = ind_info.get('max_genealo_gens', -1)
        ind_info['MinGenealGens'] = ind_info.get('min_genealo_gens', -1)
        ind_info['seed_alleles'] = ind_info.get('allele_count', -1)
        ind_info['blocks'] = ind_info.get('num_blocks', -1)
        ind_info['CentromereCount'] = ind_info.get('centromeres', -1)
        ind_info['fitness'] = ind_info.get('fitness', -1)
        ind_info['NumMuts'] = ind_info.get('mutations', -1)
        ind_info['CauseOfDeath'] = ind_info.get('cause_of_death', -1)
        info = f"{ind},{ind_info['birth_year']},{year},{ind_info['sex']},{ind_info['dad']},{ind_info['mom']},{ind_info['lifespan']},"
        info = info + f"{ind_info['lat']},{ind_info['lon']},{ind_info['marriage_state']},{ind_info['numbirths']},{ind_info['Ygens']},"
        info = info + f"{ind_info['MTGens']},{ind_info['MinGenealGens']},{ind_info['MaxGenealGens']},"
        info = info + f"{ind_info['seed_alleles']},{ind_info['CentromereCount']},{ind_info['blocks']},"
        info = info + f"{ind_info['fitness']},{ind_info['NumMuts']},{ind_info['CauseOfDeath']}"
        info = info + '\n'

    return info

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

def RIP(dead_people, IndData, model, chromosomes, mutations):

    for ind in dead_people:
        if ind in IndData:
            if IndData[ind]['marriage_state'] > -1:
                IndData[ IndData[ind]['marriage_state'] ]['marriage_state'] = -1
            del IndData[ind]
            if ind in chromosomes:
                del chromosomes[ind]
            if ind in mutations:
                del mutations[ind]

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

def initialize_population(IndData, model, free_params, chromosomes, mutations):

    if model['scenario'] == 'Eden':
        IndData = setup_pop_Eden(IndData, model, free_params, mutations)
    elif model['scenario'] == 'Flood':
        IndData = setup_pop_Flood(IndData, model, free_params, mutations)
    else:
        IndData = setup_pop_1(IndData, model, free_params, mutations)
    if model['init_heterozygosity'] > 0:
        setup_init_heterozygosity(IndData, model, free_params, chromosomes)

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

def setup_pop_Eden(IndData, model, free_params, mutations):

    for indid in range(2):
        IndData[indid] = {}
        if indid % 2 == 0:
            IndData[indid]['sex'] = 0
            IndData[indid]['marriage_state'] = indid + 1
        else:
            IndData[indid]['sex'] = 1
            IndData[indid]['marriage_state'] = indid - 1
            IndData[indid]['year_of_last_birth'] = model['spacing'] * -1
        IndData[indid]['birth_year'] = -100
        IndData[indid]['lifespan'] = 900
        IndData[indid]['fitness'] = 1
        if model['track_DNA']:
             IndData[indid]['centromeres'] = bitarray(48)
             IndData[indid]['centromeres'].setall(0)
             IndData[indid]['Y_gens'] = -1
             IndData[indid]['mt_gens'] = -1
             IndData[indid]['min_genealo_gens'] = -1
             IndData[indid]['max_genealo_gens'] = -1
             IndData[indid]['allele_count'] = 0
        if model['track_mutations']:
            numbits = free_params['numbits']
            mutations[indid] = {}
            mutations[indid][0] = [[] for _ in range(numbits)]
            mutations[indid][1] = [[] for _ in range(numbits)]
            IndData[indid]['mutations'] = 0
            IndData[indid]['fitness'] = 1
        IndData[indid]['lat'], IndData[indid]['lon'] = 0, 0

    return IndData

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

def setup_pop_Flood(IndData, model, free_params, mutations):

    for indid in range(6):
        IndData[indid] = {}
        if indid % 2 == 0:
            IndData[indid]['sex'] = 0
            IndData[indid]['marriage_state'] = indid + 1
        else:
            IndData[indid]['sex'] = 1
            IndData[indid]['marriage_state'] = indid - 1
            IndData[indid]['year_of_last_birth'] = model['spacing'] * -1
        IndData[indid]['birth_year'] = -100
        IndData[indid]['lifespan'] = 650
        IndData[indid]['fitness'] = 1
        if model['track_DNA']:
             IndData[indid]['centromeres'] = bitarray(48)
             IndData[indid]['centromeres'].setall(0)
             IndData[indid]['Y_gens'] = -1
             IndData[indid]['mt_gens'] = -1
             IndData[indid]['min_genealo_gens'] = -1
             IndData[indid]['max_genealo_gens'] = -1
             IndData[indid]['allele_count'] = 0
        if model['track_mutations']:
            numbits = free_params['numbits']
            mutations[indid] = {}
            mutations[indid][0] = [[] for _ in range(numbits)]
            mutations[indid][1] = [[] for _ in range(numbits)]
            IndData[indid]['mutations'] = 0
            IndData[indid]['fitness'] = 1
        IndData[indid]['lat'], IndData[indid]['lon'] = 0, 0

    return IndData

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

def setup_pop_1(IndData, model, free_params, mutations):

    indid = 0
    numbits = free_params['numbits']
    n = model['start_pop_size']

    structure = [0] * 100
    filename = os.path.join(data_directory, 'example_pop.csv')
    with open(filename, 'r', newline='') as file:
        reader = csv.DictReader(file)
        for row in reader:
            age = int(row['Age'])
            p = float(row['P'])
            if age < 100:
                structure[age] = int(p * n)

    for indid in range(n):
        age = random.choices(range(len(structure)), weights=structure)[0]
        IndData[indid] = {}
        IndData[indid]['sex'] = random.randint(0, 1)
        IndData[indid]['birth_year'] = -age
        IndData[indid]['lifespan'] = model['init_lifespan']
        IndData[indid]['marriage_state'] = -1
        IndData[indid]['year_of_last_birth'] = model['spacing'] * -1
        IndData[indid]['fitness'] = 1
        if model['track_DNA']:
             IndData[indid]['centromeres'] = bitarray(48)
             IndData[indid]['centromeres'].setall(0)
             IndData[indid]['Y_gens'] = -1
             IndData[indid]['mt_gens'] = -1
             IndData[indid]['min_genealo_gens'] = -1
             IndData[indid]['max_genealo_gens'] = -1
             IndData[indid]['allele_count'] = 0
        if model['track_mutations']:
            mutations[indid] = {}
            mutations[indid][0] = [[] for _ in range(numbits)]
            mutations[indid][1] = [[] for _ in range(numbits)]
            IndData[indid]['mutations'] = 0
            IndData[indid]['fitness'] = 1
        IndData[indid]['lat'], IndData[indid]['lon'] = random_coordinates()

    return IndData

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

def setup_init_heterozygosity(IndData, model, free_params, chromosomes):

    init_het = model['init_heterozygosity']
    numbits = free_params['numbits']
    for ind in IndData:
        chromosomes[ind] = {}
        chromosomes[ind][0] = bitarray(numbits)
        chromosomes[ind][1] = bitarray(numbits)
        chromosomes[ind][0].setall(0)
        chromosomes[ind][1].setall(0)
        for i in range(numbits):
            if random.random() < init_het:
                chromosomes[ind][0][i] = 1
        IndData[ind]['allele_count'] = count_alleles(ind, chromosomes)

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

def random_coordinates():

    theta = random.uniform(0, 2 * math.pi)
    r = random.uniform(0, 0.5)
    x = round(r * math.cos(theta), 2)
    y = round(r * math.sin(theta), 2)
    return x, y

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

def Innoculate_random_person(IndData, chromosomes, free_params):

    seed_individual = random.choice(list(IndData.keys()))
    free_params["seed"] = seed_individual
    if IndData[seed_individual]["sex"] == 0:
        IndData[seed_individual]["Y_gens"] = 0
    IndData[seed_individual]["mt_gens"] = 0
    IndData[seed_individual]["min_genealo_gens"] = 0
    IndData[seed_individual]["max_genealo_gens"] = 0
    chrombits = bitarray(48) # 48, because there is no chromosome 0
    chrombits.setall(1)
    IndData[seed_individual]["centromeres"] = chrombits
    numbits = free_params["numbits"]
    IndData[seed_individual]["allele_count"] = numbits * 2
    chromosomes[seed_individual] = [bitarray(numbits), bitarray(numbits)]
    chromosomes[seed_individual][0].setall(1)
    chromosomes[seed_individual][1].setall(1)

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

def count_breeding_individuals(IndData, year, model):

    count = 0
    for ind in IndData:
        age = year - IndData[ind]["birth_year"]
        sex = IndData[ind]["sex"]
        if age > model["maturity"]:
            if sex == 0:
                count += 1
            else:
                if age < IndData[ind]["lifespan"] * model["menopause"]:
                    count += 1

    return(count)

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

def calculate_misc_stats(IndData):

    Y       = sum(1 for ind in IndData.values() if ind.get("Y_gens",           -1) > 0)
    mt      = sum(1 for ind in IndData.values() if ind.get("mt_gens",          -1) > 0)
    genealo = sum(1 for ind in IndData.values() if ind.get("max_genealo_gens", -1) > 0)
    genetic = sum(1 for ind in IndData.values() if ind.get("allele_count",     -1) > 0 )
    blocks  = sum(1 for ind in IndData.values() if ind.get("num_blocks",       -1) > 0 )
    cents = sum(ind.get("centromeres", bitarray()).count() for ind in IndData.values())

    return Y, mt, genealo, genetic, blocks, cents

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

def calculate_block_stats(chromosomes):

    all_block_lengths = []
    for ind in chromosomes:
        for copy in range(2):
            bit_array = chromosomes[ind][copy]
            bit_string = bit_array.to01()
            substrings = bit_string.split("0")
            substrings = [substring for substring in substrings if substring]
            all_block_lengths.extend(map(len, substrings))
    if not all_block_lengths:
        return 0, 0, 0

    total_blocks = len(all_block_lengths)
    average_length = sum(all_block_lengths) / total_blocks
    average_length = round(average_length, 1)
    variance = round( sum((length - average_length) ** 2 for length in all_block_lengths) / total_blocks, 2)
    std_dev = variance ** 0.5
    std_dev = round(std_dev, 1)

    return total_blocks, average_length, std_dev

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

def calculate_genetic_stats(chromosomes, numbits, numinds):

    seed_genome_retained = bitarray(numbits)
    seed_genome_retained.setall(0)
    bit_counts = [0] * numbits
    tot_het = 0
    hetbits = bitarray(numbits)
    hetbits.setall(0)
    for ind in chromosomes:
        seed_genome_retained, bit_counts, tot_het = seed_counts(chromosomes, ind, numbits, seed_genome_retained, bit_counts, tot_het)
    
    numbitsretained = seed_genome_retained.count(1)
    perc_seed_genome_retained = numbitsretained / numbits * 100
    non_zero_counts = [count for count in bit_counts if count > 0]
    av_seed_genome_coverage = sum(non_zero_counts) / (len(non_zero_counts) * numbits) if non_zero_counts else 0
    av_heterozygosity = tot_het / (numbits * numinds) * 100

    return perc_seed_genome_retained, av_seed_genome_coverage, av_heterozygosity

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

def seed_counts(chromosomes, ind, numbits, seed_genome_retained, bit_counts, tot_het):

    bitarray_0 = bitarray(numbits)
    bitarray_1 = bitarray(numbits)
    bitarray_0.setall(0)
    bitarray_1.setall(0)
    if chromosomes[ind][0]:
        bitarray_0 = chromosomes[ind][0]
    if chromosomes[ind][1]:
        bitarray_1 = chromosomes[ind][1]
    seed_genome_retained |= bitarray_0 | bitarray_1
    for index in range(numbits):
        bit_counts[index] += bitarray_0[index] + bitarray_1[index]
    hetbits = bitarray_0 ^ bitarray_1
    num_het = hetbits.count(1)
    tot_het += num_het
    return seed_genome_retained, bit_counts, tot_het

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

def calculate_fitness_stats(IndData, numbits):

    numinds = len(IndData)
    NumMuts = sum(IndData[ind]['mutations'] for ind in IndData)
    AvMutsInd = NumMuts / numinds
    total_fitness = sum(IndData[ind]['fitness'] for ind in IndData)
    AvFitInd = total_fitness / numinds
    AvFitBin =  total_fitness / (numbits * numinds * 2)
    AvMutsInd = NumMuts / numinds
    AvMutsBin = NumMuts / (numbits * numinds * 2)

    return AvFitInd, AvFitBin, NumMuts, AvMutsInd, AvMutsBin

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

def load_chromosome_data(multiplier):

    filename = os.path.join(data_directory, "chromosome_data.csv")
    chromosome_arm_data = {}
    genome_size = 0
    with open(filename, 'r', newline='') as file:
        reader = csv.DictReader(file)        
        for row in reader:
            chrom = int(row['Chromosome'])
            arm_type = int(row['Arm'])         # 0 = p, 1 = q
            arm_type = 'p' if arm_type == 0 else 'q'
            arm_start = int(row['Start'] * multiplier)
            arm_length = int(row['Length'] * multiplier)
            genome_size += arm_length
            if chrom not in chromosome_arm_data:
                chromosome_arm_data[chrom] = {}
            chromosome_arm_data[chrom][arm_type] = {}
            chromosome_arm_data[chrom][arm_type]['start'] = arm_start
            chromosome_arm_data[chrom][arm_type]['length'] = arm_length
#            chromosome_arm_data[chrom][arm_type][0] = arm_start
#            chromosome_arm_data[chrom][arm_type][1] = arm_length
#    flattened_data = {chrom: dict(arm_data) for chrom, arm_data in chromosome_arm_data.items()}
#    return flattened_data, int(genome_size)
    return chromosome_arm_data, int(genome_size)

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

def load_actuarial_table():

    death_risk = {}
    filename = os.path.join(data_directory, "actuarial_table.csv")
    with open(filename, 'r', newline='') as file:
        reader = csv.DictReader(file)
        for row in reader:
            age = int(row['Age'])
            risk = float(row['Risk'])
            death_risk[age] = risk

    return death_risk

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

def setup_free_params(model):

    free_params = {}
    free_params["innoculated"] = {}
    free_params["indID"] = model["start_pop_size"] - 1
    free_params["lastpopsize"] = model["start_pop_size"]
    free_params["mutID"] = 0
    free_params["seed"] = []
    return free_params

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

def setup_output_files(model, run):

    def create_csv(filename, headers):
        with open(filename, mode='w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(headers)

    model_name = model['model_id']
    filename = os.path.join(results_directory, f"{model_name}-{run} results.csv")
    headers = ['run', 'year', 'numinds', 'marriages', 'births', 'random_deaths', 'culled_deaths', 'genetic', 'genealo', 'Y', 'mt', 'cents', 'numblocks', 'avblocksize', 'sdblocksize', 'AvFitPerInd', 'AvBinFit', 'NumMuts', 'AvMutsPerInd', 'AvMutsPerBin']
    create_csv(filename, headers)

    if model["track_dead"]:
        filename = os.path.join(results_directory, f"{model_name}-{run} deaths.csv")
        headers = ['ID', 'birthyear', 'deathyear', 'sex', 'father', 'mother', 'lifespan', 'lat', 'lon', 'married', 'numbirths', 'Ygens', 'MTgens', 'MinGenealGens', 'MaxGenealGens', 'SeedAlleles', 'CentromereCount', 'blocks', 'fitness', 'NumMuts', 'CauseOfDeath']
        create_csv(filename, headers)

    if model["mutation_hist"]:
        filename = os.path.join(results_directory, f"{model_name}-{run} mutation_histogram.csv")
        headers = ['Effect','All','EndOfRun']
        create_csv(filename, headers)

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

def setup_plot(model):

    things_to_plot_during_run = {}
    things_to_plot_during_run['year'] = []
    for variable in model["plot"]:
        things_to_plot_during_run[variable] = []
    return things_to_plot_during_run

# * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

if __name__ == "__main__":
    run_population_model()
