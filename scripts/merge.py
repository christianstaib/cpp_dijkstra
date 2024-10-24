import click
import pandas as pd


@click.command()
@click.option(
    "-p",
    "--planets_file",
    required=True,
    help="Path to the CSV file containing planets and moons data.",
)
@click.option(
    "-a",
    "--asteroids_file",
    required=True,
    help="Path to the CSV file containing asteroids data.",
)
@click.option(
    "-o",
    "--output_file",
    default="combined.csv",
    help='Output CSV file name. Defaults to "combined.csv".',
)
def combine_datasets(planets_file, asteroids_file, output_file):
    """
    Combine planets and asteroids datasets into one CSV file and adds necessary columns.
    """
    try:
        # Load planets data
        planets = pd.read_csv(planets_file, dtype=str)
        click.echo(f"Loaded planets data from {planets_file}")

        # Load asteroids data
        asteroids = pd.read_csv(asteroids_file, dtype=str)
        click.echo(f"Loaded asteroids data from {asteroids_file}")

        # Add mass and central_body columns to asteroids
        column_name = "diameter"
        col_idx = asteroids.columns.get_loc(column_name)
        asteroids.insert(col_idx + 1, "mass", "")

        column_name = "name"
        col_idx = asteroids.columns.get_loc(column_name)
        asteroids.insert(col_idx + 1, "central_body", "Sun")

        # Combine the datasets
        df = pd.concat([planets, asteroids])
        click.echo("Datasets combined successfully.")

        # Save to CSV
        df.to_csv(output_file, index=False)
        click.echo(f"Combined data saved to {output_file}")

    except Exception as e:
        click.echo(f"Error: {e}")


if __name__ == "__main__":
    combine_datasets()
