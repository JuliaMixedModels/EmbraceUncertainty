"""
Override for output of type `MixedModel`.
Without this override, a string is shown and syntax highlighting language is chosen by Highlight.jl.
With this override, the highlighting is disabled.
"""
function Books.convert_output(expr, path, out::MixedModel; kwargs...)
    return """
        ```language-plain
        $out
        ```
        """
end
