"""Integrated 250-genome scale-bearing pilot workflow.

Repeats the validated integrated SYNG -> prophage extraction/query ->
preliminary clustering/matrix workflow on the frozen 250-assembly rung.
Reuses the immutable N=100 integrated pilot (integrated-pilot-100-v1-0a11eda244a9def8)
read-only as the prior rung of record; rebuilds only release-scoped products that
legitimately depend on N. This rung is SCALE-BEARING: the ``scale_trend`` gate is
applicable (time exponent upper bound <=1.3, <=25% unexplained per-base slope change)
and drives the GO_500/NO_GO authorization for the next rung.
"""
