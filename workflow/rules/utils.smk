rule datasources_get_external:
    """Get external resources with appropriate method"""
    wildcard_constraints:
        target = wildcards_or(url_map.target)
    output:
        target = "{target}"
    input: uri = url_map.get_source
    params:
        uri = url_map.get_source_uri,
        scheme = url_map.get_source_scheme
    log: "logs/datasources_get_external/{target}.log"
    wrapper:
        f"{WRAPPER_PREFIX}/utils/datasources_get_external"


localrules: datasources_get_external
