from __future__ import annotations

import argparse
from pathlib import Path

from tract_reference_renderer.renderer import neutral_param_vector, render_tract_svg


def main() -> None:
    parser = argparse.ArgumentParser(description="Render a neutral tract SVG smoke output.")
    parser.add_argument(
        "--output",
        default="smoke_output/neutral_tract.svg",
        help="Output SVG path relative to the repo root.",
    )
    args = parser.parse_args()

    result = render_tract_svg(neutral_param_vector())
    output_path = Path(args.output)
    if not output_path.is_absolute():
        output_path = Path(__file__).resolve().parents[1] / output_path
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(result.svg, encoding="utf-8")
    print(f"wrote_svg={output_path}")
    print(f"upper_outline_points={result.diagnostics.upper_outline_points}")
    print(f"lower_outline_points={result.diagnostics.lower_outline_points}")
    print(f"tongue_outline_points={result.diagnostics.tongue_outline_points}")


if __name__ == "__main__":
    main()
