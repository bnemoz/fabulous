# API Documentation

## Overview
This document provides an overview of the API routes available in the `main.py` file. Each route is described with its endpoint, HTTP method, and a brief description of its functionality.

## Routes

### 1. `/`
- **Method:** `GET`
- **Description:** Returns a welcome message for the Fab'ulous Antibody API.

### 2. `/id`
- **Method:** `GET`, `POST`
- **Description:** Processes an antibody sequence and returns the results.
- **Parameters:**
  - `sequence_id` (string): The ID of the sequence.
  - `sequence` (string): The antibody sequence.
  - `species` (string, optional): The species of the antibody. Defaults to "human".
  - `userid` (string, optional): The user ID.
  - `authtoken` (string, optional): The authentication token.
  - `debug` (boolean, optional): If true, returns debug information.

### 3. `/ids`
- **Method:** `POST`
- **Description:** Processes multiple antibody sequences and returns the results.

### 4. `/optimize`
- **Method:** `POST`
- **Description:** Optimizes an antibody sequence.

### 5. `/clone`
- **Method:** `POST`
- **Description:** Clones an antibody sequence.

### 6. `/number`
- **Method:** `POST`
- **Description:** Numbers an antibody sequence.
- **Parameters:**
  - `header` (dict): Contains user information.
  - `data` (dict): Contains the antibody sequence and numbering scheme.

### 7. `/annotate`
- **Method:** `POST`
- **Description:** Annotates an antibody sequence.

### 8. `/clonotype`
- **Method:** `POST`
- **Description:** Identifies clonotypes in antibody sequences.

### 9. `/phylogeny`
- **Method:** `POST`
- **Description:** Generates a phylogenetic tree for antibody sequences.

### 10. `/tree`
- **Method:** `POST`
- **Description:** Builds a tree from the provided data.

### 11. `/humanize`
- **Method:** `POST`
- **Description:** Humanizes an antibody sequence.
- **Parameters:**
  - `userid` (string): The user ID.
  - `authtoken` (string): The authentication token.
  - `model` (string): The model type ("single" or "multi").
  - `debug` (boolean, optional): If true, returns debug information.
  - For `single` model:
    - `temp` (float): The temperature.
    - `sequence` (dict): Contains `sequence_id`, `sequence`, and `species`.
  - For `multi` model:
    - `oracles` (list): List of oracles.
    - `iterations` (int): Number of iterations.
    - `seqs_per_it` (int): Sequences per iteration.
    - `final_output` (string): Final output format.
    - `mutables_fwr` (int): Number of mutable framework residues.
    - `mutables_cdr` (int): Number of mutable CDR residues.
    - `initiation` (dict): Contains `sequence_id`, `sequence`, and `species`.

### 12. `/bill`
- **Method:** `POST`
- **Description:** Retrieves the billing information for a user.
- **Parameters:**
  - `userid` (string): The user ID.
  - `authtoken` (string): The authentication token.

### 13. `/graphql`
- **Method:** `POST`
- **Description:** Handles GraphQL queries.

