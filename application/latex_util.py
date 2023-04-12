import regex as re
import os

# # Для переварювання sink-ів для latex-таблиць

current_path = r'/home/goodpenguin/Desktop/loclin_shv2023/modeling/application'
os.chdir(current_path)

# print(os.getcwd())

num_examples = 5

for example_num in range(1, num_examples+1):
	current_folder = 'example' + str(example_num)
	folder_to_review = os.path.join(current_path, 'examples', current_folder)

	for filename in os.listdir(folder_to_review):
		current_file = os.path.join(folder_to_review, filename)
	
		from_file = []
	
		with open(current_file, 'r') as src:
			from_file = src.readlines()

		# print(from_file)

		pattern_for_spacings = r'(\d{2}) ([\+\- ]?\d)'
		pattern_for_start = r'^\d+ +(\d+)'

		with open(current_file, 'w') as out:
			for read_line in from_file:
				cleaned_line = re.sub(pattern_for_spacings, r'\1 & \2', read_line)
				cleaned_line = re.sub(pattern_for_start, r'\1', cleaned_line)
				
				if len(cleaned_line) > 10:
					cleaned_line = re.sub(r'\n$', r'\n\\\\\\hline\n', cleaned_line)
				
				out.write(cleaned_line)
