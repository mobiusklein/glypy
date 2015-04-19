
def reduced_end_compat(reduced_end):
    if not hasattr(reduced_end, "base_composition"):
        if len(reduced_end.links) == 0:
            reduced_end.base_composition = reduced_end.composition.clone()
        else:
            reduced_end.base_composition = reduced_end.composition.clone()
            for pos, link in reduced_end.links.items():
                reduced_end.base_composition += link.parent_loss
