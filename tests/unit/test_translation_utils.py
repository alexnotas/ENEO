import unittest

from src import translation_utils


class TestTranslationUtils(unittest.TestCase):
    def test_get_translation_missing_key_returns_fallback(self):
        value = translation_utils.get_translation("missing.key", fallback="fallback")
        self.assertEqual(value, "fallback")

    def test_set_language_invalid(self):
        current = translation_utils.get_current_language()
        result = translation_utils.set_language("xx")
        self.assertFalse(result)
        self.assertEqual(translation_utils.get_current_language(), current)

    def test_available_languages_contains_en_el(self):
        languages = translation_utils.get_available_languages()
        self.assertIn("en", languages)
        self.assertIn("el", languages)


if __name__ == "__main__":
    unittest.main()
