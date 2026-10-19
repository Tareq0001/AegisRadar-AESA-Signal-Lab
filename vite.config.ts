import { defineConfig } from 'vite';
import { resolve } from 'path';

export default defineConfig({
  base: './',
  resolve: {
    alias: {
      '@dsp': resolve(__dirname, 'src/dsp'),
      '@aesa': resolve(__dirname, 'src/aesa'),
      '@doppler': resolve(__dirname, 'src/doppler'),
    },
  },
  build: {
    target: 'es2022',
    outDir: 'dist',
    sourcemap: true,
    minify: 'esbuild',
    rollupOptions: {
      input: {
        main: resolve(__dirname, 'index.html'),
      },
    },
  },
  server: {
    port: 3000,
    open: false,
  },
});
