"""Same-set balanced union; no coordinate conversion or geometric admission."""
def balanced_union(parts, union):
    parts=tuple(parts)
    if not parts:raise ValueError('EMPTY_FIXED_PART_UNION')
    if len(parts)==1:return parts[0]
    middle=len(parts)//2
    return union(balanced_union(parts[:middle],union),balanced_union(parts[middle:],union))
