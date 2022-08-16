import multiprocessing

from collections import defaultdict

from glypy.algorithms import DistinctGlycanSet

from .graph import _enzyme_graph_inner


class Glycome(object):

    def __init__(self, glycosylases, glycosyltransferases, seeds, track_generations=False,
                 limits=None):
        if limits is None:
            limits = []
        self.glycosylases = glycosylases
        self.glycosyltransferases = glycosyltransferases
        self.seeds = seeds

        self.seen = DistinctGlycanSet()
        self.track_generations = track_generations
        self.enzyme_graph = defaultdict(_enzyme_graph_inner)
        self.history = []
        self.current_generation = DistinctGlycanSet(seeds)
        self.limits = limits

    def save_generation(self, generation):
        if self.track_generations:
            self.history.append(generation)
        self.seen.update(generation)

    def run(self, n=50):
        for i in range(n):
            generation = self.step()
            if not generation:
                break
            yield generation

    def clean_next_generation(self, generation):
        return generation - self.seen

    def within_limits(self, structure):
        for limiter in self.limits:
            if not limiter(structure):
                return False
        return True

    def step(self):
        next_generation = DistinctGlycanSet()
        for species in self.current_generation:
            parentkey = None
            for enzkey, enz in self.glycosylases.items():
                products = [root for root, leaf in enz(
                    species, refund=True) if self.within_limits(root)]
                if products:
                    if parentkey is None:
                        parentkey = str(species)
                    for product in products:
                        childkey = str(product)
                        self.enzyme_graph[parentkey][childkey].add(enzkey)
                        next_generation.add(product)
            for enzkey, enz in self.glycosyltransferases.items():
                products = [root for root in enz(
                    species) if self.within_limits(root)]
                if products:
                    if parentkey is None:
                        parentkey = str(species)
                    for product in products:
                        childkey = str(product)
                        self.enzyme_graph[parentkey][childkey].add(enzkey)
                        next_generation.add(product)
        self.save_generation(self.current_generation)
        self.current_generation = self.clean_next_generation(next_generation)
        return next_generation


def _MultiprocessingGlycome_worker(seeds_params):
    seeds, params, seen = seeds_params
    import dill
    (glycosylases, glycosyltransferases, _,
     track_generations, limits) = dill.loads(params)
    # limits.append(lambda x: x not in seen)
    glycome = Glycome(glycosylases, glycosyltransferases, seeds,
                      track_generations, limits)
    glycome.step()
    return glycome.current_generation, glycome.enzyme_graph


class MultiprocessingGlycome(Glycome):

    def __init__(self, glycosylases, glycosyltransferases, seeds, track_generations=False,
                 limits=None, processes=None):
        if processes is None:
            processes = min(multiprocessing.cpu_count(), 4)
        super(MultiprocessingGlycome, self).__init__(
            glycosylases, glycosyltransferases, seeds,
            track_generations, limits)
        self.processes = processes
        self.pool = None
        self.seen = DistinctGlycanSet()
        self._worker_params = (
            self.glycosylases, self.glycosyltransferases, tuple(),
            self.track_generations, self.limits)

    def _create_pool(self):
        self.pool = multiprocessing.Pool(self.processes)

    def _log(self, message):
        print(message)

    def log_generation_chunk(self, i, chunks, current_generation):
        self._log(".... Task %d/%d finished (%d items generated)" % (
            i, len(chunks), len(current_generation)))

    def _generate_work_loads(self, next_generation, chunks):
        import dill
        encoded_params = dill.dumps(self._worker_params)

        for i, chunk in enumerate(chunks):
            yield (chunk, encoded_params, next_generation)

    def _partition_generation(self, generation, max_chunk_size=2e3):
        n = len(generation)
        n_chunks = int(n // max_chunk_size)
        if n_chunks < self.processes * 4:
            return generation.partition(self.processes * 4)
        else:
            return generation.partition(n_chunks)

    def step(self):
        next_generation = DistinctGlycanSet()
        self._log(".... Starting Step")
        chunks = self._partition_generation(self.current_generation)
        self._log(".... Produced %d chunks" % (len(chunks) + 1,))
        work_spec = self._generate_work_loads(next_generation, chunks)
        if self.pool is None:
            self._create_pool()
        pool = self.pool
        i = 0
        for work in pool.imap_unordered(_MultiprocessingGlycome_worker, work_spec):
            i += 1
            current_generation, enzyme_graph = work
            self.log_generation_chunk(i, chunks, current_generation)
            next_generation.update(current_generation)
            for parent, children in enzyme_graph.items():
                for child, enzymes in children.items():
                    self.enzyme_graph[parent][child].update(enzymes)

        self.save_generation(self.current_generation)
        self.current_generation = self.clean_next_generation(next_generation)
        return next_generation
