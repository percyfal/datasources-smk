rule datasources_get_external:
    """Get external resources with appropriate method"""
    wildcard_constraints:
        resource_file = wildcards_or(src_map.keys())
    output:
        resource_file = "{resource_file}"
    input: uri = datasources_get_external_input
    params:
        uri = lambda wildcards: parse_uri(src_map[wildcards.resource_file])[1],
        scheme = lambda wildcards: get_uri_scheme(src_map[wildcards.resource_file])
    log: "logs/datasources_get_external/{resource_file}.log"
    wrapper:
        f"{WRAPPER_PREFIX}/utils/datasources_get_external"


localrules: datasources_get_external
