describetrees
================

Utility and "standard" for describing phylogenetic datasets (currently species tree inference only). This is like all things, heavily WIP and likely poorly designed.

See `example.toml` for an example

## Examples

Calculate average discordance (avg RF distance between true genes trees and true species trees)

```bash
describetrees example.toml --on mc1 --ref strees --against gtrees_true
```

Calculate GTEE (avg RF distance between true genes trees and true species trees)

```bash
describetrees example.toml --on mc1 --ref gtrees_true --against gtrees_est
```

Also this supports wildcards.