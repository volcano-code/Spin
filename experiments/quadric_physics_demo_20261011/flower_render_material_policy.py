"""Versioned display policy; never changes geometry or trajectory data."""

def require_material_profile(profile):
    if profile not in ('PAPER_ROLES_V1', 'MATCHED_COMPARISON_V1'):
        raise ValueError('UNKNOWN_RENDER_MATERIAL_PROFILE')
    return profile


def comparison_material_role(name, clamp_parts):
    if name in ('design', 'retained'):
        return 'polymer'
    return 'metal' if name in tuple(clamp_parts) + ('outer-race',) else 'fixed'
