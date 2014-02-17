import os
import gzip
import math

import tables as T

import numpy as NP

import plot as PL

import histogram
import utils
import nvar



def get_setting(f, name):
    """return the value of the settings with the specified name
       we use this because R replaces _ with . and this throws exceptions here
    """

    return getattr(f.root.settings, name).read()[0]





def plot_histograms(path, pair='bk', grid=None, save_to=None):


    f = T.openFile(path, mode="r")

    # figure out number of completed fparts
    fparts = f.root.histograms.ba.z.shape[0]

    # figure out number of rows and columns
    if grid is None:
        nrows = int(math.floor(math.sqrt(fparts)))
        ncols = fparts / nrows + (1 if fparts % nrows > 0 else 0)
    else:
        nrows, ncols = grid

    # select histogram
    if pair=='bk':
        hist = f.root.histograms.bk
        xlabel = '$\\beta$'
        ylabel = '$\\kappa$'
    elif pair=='ba':
        hist = f.root.histograms.ba
        xlabel = '$\\beta$'
        ylabel = '$\\alpha$'
    elif pair=='bm':
        hist = f.root.histograms.bm
        xlabel = '$\\beta$'
        ylabel = '$m$'
    elif pair=='ka':
        hist = f.root.histograms.ka
        xlabel = '$\\kappa$'
        ylabel = '$\\alpha$'
    elif pair=='km':
        hist = f.root.histograms.km
        xlabel = '$\\kappa$'
        ylabel = '$m$'
    elif pair=='am':
        hist = f.root.histograms.am
        xlabel = '$\\alpha$'
        ylabel = '$m$'

    # form histograms
    hists = [histogram.Histogram2D.initwith(hist.x, hist.y, hist.z[fpart])
             for fpart in xrange(fparts)]

    # form titles
    steps_per_fpart = f.root.settings.steps[0] / f.root.settings.fparts[0]
    titles = ['step=%d' % ((fpart+1) * steps_per_fpart) for fpart in xrange(fparts-1)]
    titles.append('step=%d' % len(f.root.dynamics.intra.M))

    # form xticks and yticks
    ## max_beta = get_setting(f, 'max.beta')
    ## max_alpha = get_setting(f, 'max.alpha')

    ## xticks = [("0", 0), ("%.0f" % max_beta, max_beta)]
    ## yticks = [("0", 0), ("%.0f" % max_alpha, max_alpha)]

    fig = histogram.Histogram2D.plot_many(hists, nrows, ncols,
                                          #xticks=xticks, yticks=yticks,
                                          xlabel=xlabel, ylabel=ylabel,
                                          ##xlabel='$\\beta$', ylabel='$\\alpha$',
                                          titles=titles,
                                          grid=True, titlesize='large')
    
    f.close()

    if save_to is None:
        fig.show()
        del(fig)
    else:
        fig.save(save_to)
    






#===========================================================================
#                          POPULATION DESIGN
#===========================================================================


class Profile:

    pid = 0

    def __init__(self, beta, kappa, alpha):

        self.beta = beta
        self.kappa = kappa
        self.alpha = alpha
        self.cn = 0
        self.pid = Profile.pid
        Profile.pid += 1

    def inc(self, step):
        self.cn += step

    def dec(self, step):
        self.cn -= step
        
    def __str__(self):
        return 'b=%.3f|k=%.3f|a=%.3f' % (self.beta, self.kappa, self.alpha)

    def tostring(self):
        return "(%d 0 %.3f %.3f %.3f %d 0 0)" % \
               (self.pid, self.beta, self.kappa, self.alpha, self.cn)



class Plasmid:

    def __init__(self, profile, cn):

        self.profile = profile
        self.cn = cn

    def tostring(self):
        return "(%d %d)" % (self.profile.pid, self.cn)



class Host:

    def __init__(self, omega_0=1):
        self.omega = 1 + NP.random.uniform()
        self.pool = {} # a hash with keys:Profile, values:Plasmid

    def add_plasmid(self, plasmid):
        try:
            self.pool[plasmid.profile].cn += plasmid.cn
        except KeyError:
            # profile does not exist in the pool -- add a new entry
            if plasmid.cn > 0:
                self.pool[plasmid.profile] = plasmid

    def tostring(self):
        return "(%.3e 0.05 1 0 0 (%s))" % \
               (self.omega, ' '.join([p.tostring() for p in self.pool.values()]))




def construct_population(fname_out, wt_tuple, mut_tuple, n=1000,
                         r_wt=7, r_mut=0.01, print_info=False):
    """
    construct a population and save it in FNAME_OUT (and HDF file)
    (WT|MUT)_TUPLE are 3-tuples containing (b,k,a) profile configurations
    N specifies the population size
    R_WT is the average number of the wild-type plasmids per host
    R_MUT is the fraction of mutants with respect to the WT
    """

    # reset Profile pid
    Profile.pid = 0

    # construct profiles
    wt_profile = Profile(*wt_tuple)
    mut_profile = Profile(*mut_tuple)

    # determine number of WT plasmids in each host
    n_wt = NP.random.poisson(r_wt, n)
    
    # calculate number of MUTANT plasmids
    # based on the number of WT ones
    n_mut = int(NP.ceil(r_mut * NP.sum(n_wt)))

    if (print_info):
        print 'total WT=', NP.sum(n_wt)
        print 'MUT=', n_mut

    # assemble hosts
    hosts = []

    for i in xrange(n):
        host = Host()
        host.add_plasmid(Plasmid(wt_profile, n_wt[i]))
        hosts.append(host)

    # distribute mutants among hosts
    for i in xrange(n_mut):
        # choose a random host and add a single mutant copy
        hosts[NP.random.randint(len(hosts))].add_plasmid(
            Plasmid(mut_profile, 1))

    # form host string
    hosts_st = '(' + ' '.join([h.tostring() for h in hosts]) + ')'
    # form profile string
    prof_st = '(%s %s)' % (wt_profile.tostring(), mut_profile.tostring())

    # write to HDF file
    fout = T.openFile(fname_out, 'w')

    # create population group
    pgroup = fout.createGroup("/", "population")

    # create profile and hosts string
    fout.createArray(pgroup, "profiles", prof_st)
    fout.createArray(pgroup, "hosts", hosts_st)
    
    fout.close()






# BETA-ALPHA EXPERIMENTS
def generate_populations_ba(path="populations/ba",
                            wt_tuple=(0.4, 0.9, 0.9),
                            beta_values=NP.arange(0.35, 0.51, 0.01),
                            alpha_values=NP.arange(0.85, 1.01, 0.01),
                            r_wt=7, r_mut=0.01):

    # if PATH does not exist -- create it
    if not os.path.exists(path):
        os.makedirs(path)

    for i, beta in enumerate(beta_values):
        for j, alpha in enumerate(alpha_values):
            fname = os.path.join(path, "population.%d.%d.h5" % (i,j))
            construct_population(fname, wt_tuple,
                                 (beta, wt_tuple[1], alpha),
                                 r_wt=r_wt, r_mut=r_mut)

    print "BA: Created %d populations in %s" % (len(beta_values)*len(alpha_values), path)





# KAPPA-ALPHA EXPERIMENTS
def generate_populations_ka(path="populations/ka",
                            wt_tuple=(0.4, 0.9, 0.9),
                            kappa_values=NP.arange(0.80, 0.951, 0.01),
                            alpha_values=NP.arange(0.85, 1.0, 0.01),
                            r_wt=7, r_mut=0.01):

    # if PATH does not exist -- create it
    if not os.path.exists(path):
        os.makedirs(path)

    for i, kappa in enumerate(kappa_values):
        for j, alpha in enumerate(alpha_values):
            fname = os.path.join(path, "population.%d.%d.h5" % (i,j))
            construct_population(fname, wt_tuple,
                                 (wt_tuple[0], kappa, alpha),
                                 r_wt=r_wt, r_mut=r_mut)


    print "KA: Created %d populations in %s" % (len(kappa_values)*len(alpha_values), path)



#                          ============  
#                          === TODO ===
#                          ============  

# === process the COMP results and draw a contour plot with the results ===


# ===================================================================
#                        AD-HOC FUNCTIONS
# ===================================================================


def rename_files(path, doit=False):
    "convert results.%d.%d.%d.h5 to results.%d.0.%d.%d.h5"

    fnames = os.listdir(path)

    for fname in fnames:

        tokens = fname.split('.')

        if (len(tokens) == 5 and tokens[0] == "results" and tokens[-1] == "h5"):
            new_fname = '.'.join(tokens[:2]) + '.0.' + '.'.join(tokens[2:])
            if (doit):
                os.rename(os.path.join(path, fname),
                          os.path.join(path, new_fname))
            print '%s ==> %s' % (fname, new_fname)




def calc_fitness(cn, cn_opt, sigma, alpha, gamma):
    "calc and return the fitness of a host"

    #return NP.sqrt(1-alpha) * NP.exp(- ((cn - cn_opt) / sigma) ** 2 / 2.)

    #return NP.exp(- ((cn - cn_opt) / sigma) ** 2 / 2.)

    #return (1 - NP.exp(-0.5*cn*alpha)) * NP.exp(- ((cn - cn_opt) / sigma) ** 2 / 2.)

    #return (1 - cn*alpha / (10. + cn*alpha)) * NP.exp(- ((cn - cn_opt) / sigma) ** 2 / 2.)

    return NP.exp( - ( (cn - cn_opt) / sigma) ** 2 / 2.) - gamma * cn * alpha

    #return NP.exp(- ((cn_opt - (1+alpha) * cn) / sigma) ** 2 / 2.)

    #return NP.exp(- ((cn_opt - (1+alpha) * cn) / sigma) ** 2 / 2.) / (sigma * NP.sqrt(2*NP.pi))





def plot_fitness(cn=10, sigma=5, gamma=0.1):

    n = NP.arange(0, 40, 0.01)

    print 'f(0) = %.3f' % calc_fitness(0, cn, sigma, 0, 0)

    #A = [0, 0.5, 1, 1.5]
    A = [0, 0.1, 0.5]

    fig = PL.Figure()

    fig.plot(n,[calc_fitness(n, cn, sigma, a, gamma) for a in A],
             vlines=[(cn, None, 'k:')], hlines=[(None, 0, 'k-')],
             #vlines=[(cn - gamma*a*sigma**2, None, 'k:') for a in A],
             legends=['a=%.1f' % a for a in A])

    
    slopes = [- (n - cn) * NP.exp(- ((n-cn) / sigma) ** 2 / 2.) / (sigma ** 2) - gamma * alpha
              for alpha in A]
    fig.plot(n, slopes)

    fig.show()





    

def plot_hill(cn, psi=1):

    m = NP.arange(0, 1, 0.01)

    PL.Figure().plot(m, ((cn*m) / (psi + cn*m), (cn*m) / (10*psi + cn*m)),
                     xlabel='m', title='$\\sum nm / (\\psi + \\sum nm)$').show()

