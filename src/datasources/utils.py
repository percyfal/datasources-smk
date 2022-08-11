##################################################
# Formatting functions and other utilities
##################################################
def wildcards_or(items):
    items = [str(x) for x in set(items)]
    return f'({"|".join(items)})'
